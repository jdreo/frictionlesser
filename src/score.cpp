#include <cmath>
#include <vector>

#include "frictionless/frictionless.h"
#include "frictionless/score.h"
#include "R/pgamma.h"

namespace frictionless {

FriedmanScore::FriedmanScore( const Transcriptome& rt, const double a) :
    _transcriptome(rt),
    alpha(a),
    _has_init_signature(false),
    _has_size_cache(false),
    _has_transcriptome_cache(false)
{
    ASSERT(_transcriptome.ranks( ).size() > 0);
    ASSERT(_transcriptome.ranks(0).size() > 0);
    ASSERT(_transcriptome.genes_nb() > 0);
    ASSERT(_transcriptome.samples_nb() > 0);
    ASSERT(_transcriptome.cells_nb(0) > 0);
    ASSERT(alpha > 0);

    const size_t samples_nb = _transcriptome.samples_nb();

    _transcriptome_cache.reserve(samples_nb); // E, F, GG, T, SSR
    _size_cache.reserve(samples_nb); // B, C
    _size_cache.signature_size = 0;
    _swap_cache.reserve(samples_nb, _transcriptome.cells_nb()); // R, A, D
}

void FriedmanScore::clear_cache()
{
    CLUTCHLOG(debug, "Clear cache");
    _transcriptome_cache.clear(); // E, F, GG, T, SSR
    _has_transcriptome_cache = false;
    _swap_cache.clear(); // R, A, D
    _size_cache.clear();
    _size_cache.signature_size = 0;
    _has_size_cache = false;

    _cached_gene_in = 0;
    _cached_gene_out = 0;
}

/******************************************
 * E, F, GG, SSR, T
 ******************************************/
void FriedmanScore::new_transcriptome(const bool print_progress)
{
    CLUTCHLOG(progress, "New transcriptome...");
    clear_cache();

    CLUTCHLOG(info,"Precompute constants");
    double progress = 0;
    for(size_t i : _transcriptome.samples() ) { // All samples.
        if(print_progress) {
            std::clog << "\r" << static_cast<size_t>(progress / _transcriptome.samples_nb() * 100.0) << "%     ";
            std::clog.flush();
            CLUTCHCODE(debug,
                std::clog << std::endl;
           );
            progress++;
        }
        CLUTCHLOGD(debug, "Sample: " << i, 1);
        // Constants depending only on the number of cells
        // in each sample: _transcriptome_cache.E, _transcriptome_cache.F, _transcriptome_cache.GG.
        const double m_i = _transcriptome.cells_nb(i);
        _transcriptome_cache.E .push_back( 3 * m_i * std::pow(m_i+1,2) );
        _transcriptome_cache.F .push_back( m_i * (m_i+1) );
        _transcriptome_cache.GG.push_back( m_i / (m_i - 1) );
        CLUTCHLOGD(xdebug, "E=" << _transcriptome_cache.E.back() << ",    F=" << _transcriptome_cache.F.back() << ",    GG=" << _transcriptome_cache.GG.back(), 2);

        // Constants for each cells of each sample: _transcriptome_cache.SSR, _transcriptome_cache.T.
        std::vector<double> SSR_j;
        std::vector<double> T_j;
        for(size_t j : _transcriptome.genes()) { // All geneset.

            double sum_ranks_sq = 0; // for _transcriptome_cache.SSR.
            // Rank => ties count^3, for _transcriptome_cache.T.
            // This should work even if ranks are unordered.
            std::map<double,unsigned long long> ties; // For _transcriptome_cache.T.
            for(size_t c : _transcriptome.cells(i) ) { // All cells of this sample.
                const double rk = _transcriptome.rank(c,j);
                sum_ranks_sq += std::pow(rk,2); // For _transcriptome_cache.SSR.
                if(ties.contains(rk)) { // For _transcriptome_cache.T.
                    ties[rk]++;
                } else {
                    ties[rk] = 1;
                }
                ASSERT(ties[rk]!=0); // Overflow guard.
            } // for c in cells
            SSR_j.push_back( sum_ranks_sq );

            // Sum cubic number of ties.
            double sum_ties_cube = 0;
            for([[maybe_unused]] const auto& [r,t] : ties) {
                sum_ties_cube += std::pow(t,3);
            }
            T_j.push_back( sum_ties_cube );

        } // for j in geneset
        _transcriptome_cache.SSR.push_back(SSR_j);
        _transcriptome_cache.T.push_back(T_j);
        CLUTCHCODE(xdebug,
            std::ostringstream ssr;
            for(double s : SSR_j) { ssr << " " << s; }
            CLUTCHLOGD(xdebug, "SSR={" << ssr.str() << "}", 2);
            std::ostringstream ts;
            for(double t : T_j  ) {  ts << " " << t; }
            CLUTCHLOGD(xdebug, "T={" << ts.str() << "}", 2);
        );

    } // for i in samples
    if(print_progress) {std::clog << std::endl;}

#ifndef NDEBUG
    // Basic size checks.
    ASSERT(_transcriptome_cache.E  .size() == _transcriptome.samples_nb());
    ASSERT(_transcriptome_cache.F  .size() == _transcriptome.samples_nb());
    ASSERT(_transcriptome_cache.GG .size() == _transcriptome.samples_nb());
    ASSERT(_transcriptome_cache.SSR.size() == _transcriptome.samples_nb());
    ASSERT(_transcriptome_cache.T  .size() == _transcriptome.samples_nb());
    for(auto row : _transcriptome_cache.SSR) {
        ASSERT(row.size() == _transcriptome.genes_nb()); }
    for(auto row : _transcriptome_cache.T) {
        ASSERT(row.size() == _transcriptome.genes_nb()); }
#endif
    _has_transcriptome_cache = true;
    CLUTCHLOG(note, "OK");
}

/******************************************
 * B, C
 ******************************************/
void FriedmanScore::new_signature_size(const size_t signature_size)
{
    CLUTCHLOG(progress, "New signature size: " << signature_size << "...");
    ASSERT(_has_transcriptome_cache);
    CLUTCHLOG(xdebug, "Cached signature size: " << _size_cache.signature_size);
    _size_cache.clear();
    ASSERT(signature_size > 0);
    _size_cache.signature_size = signature_size;

    for(size_t i : _transcriptome.samples()) {
        ASSERT(_transcriptome.cells_nb(i) > 0);
        #ifndef NDEBUG
            _size_cache.B.push_back( std::pow(signature_size,2) * _transcriptome_cache.E[i] );
            _size_cache.C.push_back(          signature_size    * _transcriptome_cache.F[i]);
        #else
            _size_cache.B.push_back( std::pow(signature_size,2) * _transcriptome_cache.E.at(i) );
            _size_cache.C.push_back(          signature_size    * _transcriptome_cache.F.at(i));
        #endif
        CLUTCHLOGD(xdebug, "Sample: " << i, 1);
        CLUTCHLOGD(xdebug, "_size_cache.B=" << _size_cache.B.back() << ",    _size_cache.C=" << _size_cache.C.back(), 2);
    } // for i in samples

    _has_size_cache = true;
    CLUTCHLOG(note, "OK");
}

/******************************************
 * R, A, D
 ******************************************/
// void FriedmanScore::new_swap(const size_t gene_in, const size_t gene_out)
// {
//     CLUTCHLOG(debug, "New swap: in=" << gene_in << ", out=" << gene_out);
//     ASSERT(gene_in != gene_out);
//     ASSERT(not (gene_in == _cached_gene_in and gene_out == _cached_gene_out));

//     for(size_t i : _transcriptome.samples()) {
//         double sum_Rin = 0;
//         double sum_r2 = 0;
//         for(size_t c : _transcriptome.cells(i)) {
//             // FIXME double check that
//             // R_c(G\setminus j_{in}) == R_c(G)-r_{c,j_{in}}
//             #ifndef NDEBUG
//                 sum_Rin +=  _transcriptome.rank(c,gene_out) * (_swap_cache.R[c] - _transcriptome.rank(c,gene_in))
//                           - _transcriptome.rank(c,gene_in ) *  _swap_cache.R[c];
//             #else
//                 sum_Rin +=  _transcriptome.rank(c,gene_out) * (_swap_cache.R.at(c) - _transcriptome.rank(c,gene_in))
//                           - _transcriptome.rank(c,gene_in ) *  _swap_cache.R.at(c);
//             #endif
//             sum_r2 += std::pow(_transcriptome.rank(c,gene_in),2)
//                     + std::pow(_transcriptome.rank(c,gene_out),2);
//         } // for c in cells
//         #ifndef NDEBUG
//             const double Ai = _swap_cache.A[i];
//             _swap_cache.A[i] = Ai + 24*sum_Rin + 12*sum_r2;

//             const double Di = _swap_cache.D[i];
//             _swap_cache.D[i] = Di + 1/(_transcriptome.cells_nb(i)-1)
//                  * ( _transcriptome_cache.T[i][gene_out] - _transcriptome_cache.T[i][gene_in] );
//          #else
//             const double Ai = _swap_cache.A.at(i);
//             _swap_cache.A.at(i) = Ai + 24*sum_Rin + 12*sum_r2;

//             const double Di = _swap_cache.D.at(i);
//             _swap_cache.D.at(i) = Di + 1/(_transcriptome.cells_nb(i)-1)
//                  * ( _transcriptome_cache.T.at(i).at(gene_out) - _transcriptome_cache.T.at(i).at(gene_in) );
//          #endif
//     } // for i in samples

//     _cached_gene_in = gene_in;
//     _cached_gene_out = gene_out;
// }
void FriedmanScore::new_swap(const size_t gene_in, const size_t gene_out)
{
    CLUTCHLOG(xdebug, "New swap: in=" << gene_in << ", out=" << gene_out);
    ASSERT(_has_transcriptome_cache);
    ASSERT(_has_size_cache);
    ASSERT(_has_init_signature);
    ASSERT(gene_in != gene_out);
    ASSERT(not (gene_in == _cached_gene_in and gene_out == _cached_gene_out));

    for(size_t i : _transcriptome.samples()) {
        for(size_t c : _transcriptome.cells(i)) {
            // Substract gene_out.
            _swap_cache.A[i] -= 24 * _swap_cache.R[c] * _transcriptome.rank(c,gene_out);
            _swap_cache.R[c] -= _transcriptome.rank(c,gene_out);
            // Add gene_in.
            _swap_cache.A[i] += 24 * _swap_cache.R[c] * _transcriptome.rank(c,gene_in);
            _swap_cache.R[c] += _transcriptome.rank(c,gene_in);
        } // for c in cells

        _swap_cache.A[i] += 12 * _transcriptome_cache.SSR[i][gene_in];
        _swap_cache.A[i] += 12 * _transcriptome_cache.SSR[i][gene_out];

        const size_t m_i = _transcriptome.cells_nb(i);
        _swap_cache.D[i] += _transcriptome_cache.T[i][gene_in]  / (m_i - 1);
        _swap_cache.D[i] -= _transcriptome_cache.T[i][gene_out] / (m_i - 1);
    } // for i in samples

    _cached_gene_in = gene_in;
    _cached_gene_out = gene_out;
}

/******************************************
 * R, A, D
 ******************************************/
void FriedmanScore::init_signature(const Signature& geneset)
{
    CLUTCHLOG(debug, "Initialize signature: " << geneset.str());
    ASSERT(_has_transcriptome_cache);
    ASSERT(_has_size_cache);
    ASSERT(geneset.selected.size()+geneset.rejected.size() == _transcriptome.genes_nb());
    ASSERT(_size_cache.signature_size == geneset.selected.size());

    _swap_cache.clear(); // R, A, D
    CLUTCHLOGD(xdebug, "R:", 1);
    for(size_t c : _transcriptome.cells()) {
        double sum_r = 0;
        for(size_t j : geneset.selected) {
            sum_r += _transcriptome.rank(c,j);
        }
        _swap_cache.R.push_back(sum_r);
        CLUTCHLOGD(xdebug, "R_" << c << "=" << _swap_cache.R.back(), 2);
    } // for c in cells

    for(size_t i : _transcriptome.samples()) {
        CLUTCHLOGD(xdebug, "Sample: " << i, 1);
        double sum_Rc2 = 0;
        for(size_t c : _transcriptome.cells(i)) {
            #ifndef NDEBUG
                sum_Rc2 += 12*std::pow(_swap_cache.R[c], 2);
            #else
                sum_Rc2 += 12*std::pow(_swap_cache.R.at(c), 2);
            #endif
        } // c in sample cells
        // _swap_cache.A.push_back(12 * sum_Rc2);
        _swap_cache.A.push_back(sum_Rc2);

        // Paper version:
        // double sum_t = 0;
        // const double m_i = _transcriptome.cells_nb(i);
        // for(size_t j : _transcriptome.geneset()) {
        //     #ifndef NDEBUG
        //         const bool selected = geneset[j];
        //         sum_t += selected * (_transcriptome_cache.T[i][j] - m_i);
        //     #else
        //         const bool selected = geneset.at(j);
        //         sum_t += selected * (_transcriptome_cache.T.at(i).at(j) - m_i);
        //     #endif
        // }
        // _swap_cache.D.push_back( 1/(m_i-1) * sum_t );

        // Zakiev code:
        double sum_t = 0;
        for(size_t j : geneset.selected) {
            sum_t += _transcriptome_cache.T[i][j];
        }
        const double m_i = _transcriptome.cells_nb(i);
        _swap_cache.D.push_back( (1/(m_i-1)) * sum_t - _transcriptome_cache.GG[i] * _size_cache.signature_size );

        CLUTCHLOGD(xdebug, "A=" << _swap_cache.A.back() << ",    D=" << _swap_cache.D.back(), 2);

    } // for i in samples

    _has_init_signature = true;
}


double FriedmanScore::score(const Signature& geneset)
{
    CLUTCHLOG(xdebug, "Compute score of: " << geneset);
    ASSERT(_has_transcriptome_cache);
    ASSERT(_has_size_cache);
    ASSERT(_has_init_signature);

    // const size_t current_signature_size = geneset.selected.size();
    // ASSERT(current_signature_size > 0);
    // if(_size_cache.signature_size != current_signature_size) {
    //     CLUTCHLOG(warning, "Current signature size " << current_signature_size
    //         << " does not match cached signature size of " << _size_cache.signature_size
    //         << ", I will silently call for a new_signature_size");
    //     new_signature_size(current_signature_size);
    // }

    double score = 0;
    for(size_t i : _transcriptome.samples()) {
        CLUTCHLOGD(xdebug, "Sample: " << i, 1);
        double s_hat;
        #ifndef NDEBUG
            const double CD = _size_cache.C[i] - _swap_cache.D[i];
        #else
            const double CD = _size_cache.C.at(i) - _swap_cache.D.at(i);
        #endif
        if(CD != 0) {
            #ifndef NDEBUG
                s_hat = ( (_swap_cache.A[i]-_size_cache.B[i]) / CD )
            #else
                s_hat = ( (_swap_cache.A.at(i)-_size_cache.B.at(i)) / CD )
            #endif
                  / std::pow(geneset.selected.size(), this->alpha);
        } else {
            s_hat = 0;
        }
        CLUTCHLOGD(xdebug, "s_hat=" << s_hat, 2);
        ASSERT(not std::isnan(s_hat));
        ASSERT(s_hat >= 0);
        ASSERT(not std::isinf(s_hat));

        score += sqrt_logchisq(s_hat, _transcriptome.cells_nb(i)-1);
    }

    CLUTCHLOGD(xdebug, "Score: " << score, 1);
    ASSERT(not std::isnan(score));
    ASSERT(score >= 0);
    ASSERT(not std::isinf(score));
    return score;
}

double FriedmanScore::sqrt_logchisq(const double s, const double m) const
{
    return std::sqrt(-1 * R::pgamma(/*x*/s, /*alph*/m/2, /*scale*/2, /*lower_tail*/0, /*log_p*/1));
}

bool FriedmanScore::has_init_signature()
{
    return _has_init_signature;
}

bool FriedmanScore::has_transcriptome_cache()
{
    return _has_transcriptome_cache;
}

bool FriedmanScore::has_size_cache()
{
    return _has_size_cache;
}

CacheSwap& FriedmanScore::swap_cache()
{
    return _swap_cache;
}

void FriedmanScore::swap_cache( const CacheSwap& cache)
{
    _swap_cache = cache;
}

void FriedmanScore::swap_cache( CacheSwap&& cache)
{
    _swap_cache = std::move(cache);
}

void FriedmanScore::load_transcriptome_cache(std::istream& in)
{
    this->_transcriptome_cache.load(in);

    // Basic size checks.
    ASSERT(_transcriptome_cache.E  .size() == _transcriptome.samples_nb());
    ASSERT(_transcriptome_cache.F  .size() == _transcriptome.samples_nb());
    ASSERT(_transcriptome_cache.GG .size() == _transcriptome.samples_nb());
    ASSERT(_transcriptome_cache.SSR.size() == _transcriptome.samples_nb());
    ASSERT(_transcriptome_cache.T  .size() == _transcriptome.samples_nb());
    for(auto row : _transcriptome_cache.SSR) {
        ASSERT(row.size() == _transcriptome.genes_nb()); }
    for(auto row : _transcriptome_cache.T) {
        ASSERT(row.size() == _transcriptome.genes_nb()); }

    _has_transcriptome_cache = true;
}


void FriedmanScore::load_size_cache(std::istream& in)
{
    this->_size_cache.load(in);

    // Basic size checks.
    ASSERT(_size_cache.B.size() == _transcriptome.samples_nb());
    ASSERT(_size_cache.C.size() == _transcriptome.samples_nb());
    ASSERT(_size_cache.signature_size > 0);

    _has_size_cache = true;
}

void FriedmanScore::save_transcriptome_cache(std::ostream& out)
{
    // Basic size checks.
    ASSERT(_transcriptome.samples_nb() > 0);
    ASSERT(_transcriptome.genes_nb() > 0);
    ASSERT(_transcriptome_cache.E  .size() == _transcriptome.samples_nb());
    ASSERT(_transcriptome_cache.F  .size() == _transcriptome.samples_nb());
    ASSERT(_transcriptome_cache.GG .size() == _transcriptome.samples_nb());
    ASSERT(_transcriptome_cache.SSR.size() == _transcriptome.samples_nb());
    ASSERT(_transcriptome_cache.T  .size() == _transcriptome.samples_nb());
    for(auto row : _transcriptome_cache.SSR) {
        ASSERT(row.size() == _transcriptome.genes_nb()); }
    for(auto row : _transcriptome_cache.T) {
        ASSERT(row.size() == _transcriptome.genes_nb()); }

    this->_transcriptome_cache.save(out);
}

void FriedmanScore::save_size_cache(std::ostream& out)
{
    // Basic size checks.
    ASSERT(_transcriptome.samples_nb() > 0);
    ASSERT(_size_cache.B.size() == _transcriptome.samples_nb());
    ASSERT(_size_cache.C.size() == _transcriptome.samples_nb());
    ASSERT(_size_cache.signature_size > 0);

    this->_size_cache.save(out);
}

} // frictionless
