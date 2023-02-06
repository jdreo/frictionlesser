#include <frictionless/frictionless.h>
#include <frictionless/cache.h>

namespace frictionless {

void Serialize::save(std::ostream& out, const std::vector<double>& vec, const size_t depth) const
{
    CLUTCHLOGD(xdebug, "save vector of size: " << vec.size(), depth);
    this->save(out, vec.size(), depth+1);
    for(const double& v : vec) {
        this->save(out, v, depth+1);
    }
}

void Serialize::save(std::ostream& out, const std::vector<std::vector<double>>& table, const size_t depth) const
{
    CLUTCHLOGD(xdebug, "save table of size: " << table.size(), depth);
    this->save(out, table.size(), depth+1);
    for(const std::vector<double>& vec : table) {
        this->save(out, vec, depth+1);
    }
}

std::vector<double> Serialize::load_vector(std::istream& in, const size_t depth) const
{
    CLUTCHLOGD(xdebug, "load vector", depth);
    auto size = this->load_scalar<std::vector<double>::size_type>(in);

    std::vector<double> vec;
    vec.reserve(size);
    for(size_t i = 0; i < size; ++i) {
        vec.push_back( this->load_scalar<double>(in, depth+1) );
    }
    return vec;
}

std::vector<std::vector<double>> Serialize::load_table(std::istream& in, const size_t depth) const
{
    CLUTCHLOGD(xdebug, "load table", depth);
    auto size = this->load_scalar<std::vector<std::vector<double>>::size_type>(in);

    std::vector<std::vector<double>> table;
    table.reserve(size);
    for(size_t i = 0; i < size; ++i) {
        table.push_back( this->load_vector(in, depth+1) );
    }
    return table;
}


void CacheTranscriptome::save(std::ostream& out) const
{
    CLUTCHLOG(info, "Save transcriptome cache...");
    ASSERT(E  .size() > 0);
    ASSERT(F  .size() > 0);
    ASSERT(GG .size() > 0);
    ASSERT(T  .size() > 0);
    ASSERT(SSR.size() > 0);
    ASSERT(E  .size() == F  .size());
    ASSERT(F  .size() == GG .size());
    ASSERT(GG .size() == T  .size());
    ASSERT(T  .size() == SSR.size());

    Serialize::save(out, E);
    Serialize::save(out, F);
    Serialize::save(out, GG);

    Serialize::save(out, T);
    Serialize::save(out, SSR);
    CLUTCHLOG(info, "OK");
}

void CacheTranscriptome::load(std::istream& in)
{
    CLUTCHLOG(info, "Load transcriptome cache...");
    clear();

    E   = Serialize::load_vector(in);
    F   = Serialize::load_vector(in);
    GG  = Serialize::load_vector(in);

    T   = Serialize::load_table(in);
    SSR = Serialize::load_table(in);

    ASSERT(E  .size() > 0);
    ASSERT(F  .size() > 0);
    ASSERT(GG .size() > 0);
    ASSERT(T  .size() > 0);
    ASSERT(SSR.size() > 0);
    ASSERT(E  .size() == F  .size());
    ASSERT(F  .size() == GG .size());
    ASSERT(GG .size() == T  .size());
    ASSERT(T  .size() == SSR.size());
    CLUTCHLOG(info, "OK");
}

void CacheTranscriptome::reserve(const size_t samples_nb)
{
    E  .reserve(samples_nb);
    F  .reserve(samples_nb);
    GG .reserve(samples_nb);
    T  .reserve(samples_nb);
    SSR.reserve(samples_nb);
}

void CacheTranscriptome::clear()
{
    CLUTCHLOG(debug, "Clear transcriptome cache");
    E  .clear();
    F  .clear();
    GG .clear();
    T  .clear();
    SSR.clear();
}


void CacheSize::save(std::ostream& out) const
{
    CLUTCHLOG(info, "Save size cache...");
    ASSERT(signature_size > 0);
    ASSERT(B.size() > 0);
    ASSERT(C.size() > 0);
    ASSERT(B.size() == C.size());

    Serialize::save(out, signature_size);
    Serialize::save(out, B);
    Serialize::save(out, C);
    CLUTCHLOG(info, "OK");
}

void CacheSize::load(std::istream& in)
{
    CLUTCHLOG(info, "Load size cache...");
    clear();

    signature_size = Serialize::load_scalar<size_t>(in);
    B = Serialize::load_vector(in);
    C = Serialize::load_vector(in);

    ASSERT(signature_size > 0);
    ASSERT(B.size() > 0);
    ASSERT(C.size() > 0);
    ASSERT(B.size() == C.size());
    CLUTCHLOG(info, "OK");
}

void CacheSize::reserve(const size_t samples_nb)
{
    B.reserve(samples_nb);
    C.reserve(samples_nb);
}

void CacheSize::clear()
{
    CLUTCHLOG(info, "Clear size cache.");
    signature_size = 0;
    B.clear();
    C.clear();
}


CacheSwap::CacheSwap( // copy
    const std::vector<double>& r,
    const std::vector<double>& a,
    const std::vector<double>& d
    ) : R(r), A(a), D(d), _has_cache(true)
{}

//! Move constructor.
CacheSwap::CacheSwap( // move
    std::vector<double>&& r,
    std::vector<double>&& a,
    std::vector<double>&& d
    ) : R(std::move(r)), A(std::move(a)), D(std::move(d)), _has_cache(true)
{}

CacheSwap::CacheSwap() : _has_cache(false) {}

//! Reserve memory for all members.
void CacheSwap::reserve(const size_t samples_nb, const size_t cells_nb)
{
    A.reserve(samples_nb);
    D.reserve(samples_nb);
    R.reserve(cells_nb);
}

//! Empty all cache.
void CacheSwap::clear()
{
    R.clear();
    A.clear();
    D.clear();
}



} // frictionless
