

#ifndef CONFIG_H
#define CONFIG_H

#include <cmath>
#include <assert.h>
#include <iostream>
#include <vector>

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

struct Pose {
    std::vector<double> v;
    Pose () {}
};

class Config {
public:
    Config ();
    Config (const Config &c);
    Config (int dim_t, int dim_r, int source=-1, bool ws=true);
    Config (int dim_t, std::vector<double> t, int dim_r, std::vector<double> r, int source=-1, bool ws=true);
    ~Config ();

    void init (int tt, int rr, bool ws=true) {
        dim_t = tt;
        dim_r = rr;
        t.resize(dim_t);
        r.resize(dim_r);
        this->ws = ws;
    }
    void init (const int &src_dim, const std::vector<double> &src, int &tgt_dim, std::vector<double> &tgt) {
        tgt_dim = src_dim;
        tgt.resize(src_dim);
    }

    void copy (const int &src_dim, const std::vector<double> &src, int &tgt_dim, std::vector<double> &tgt) {
        init(src_dim, src, tgt_dim, tgt);
        for (int i=0;i<tgt_dim;++i)
            tgt[i] = src[i];
    }

    inline double normT () const {
        //assert(ws);
        double nt(0);
        for (int i=0;i<dim_t;++i)
            nt += SQR(t[i]);
        return sqrt(nt);
    }

    inline double normR () const {
        //assert(ws);
        double nr(0);
        for (int i=0;i<dim_r;++i)
            nr += SQR(r[i]);
        return sqrt(nr);
    }

    inline double norm () const {
        assert(ws);
        return sqrt(normsqr());
    }

    inline double normsqr () const {
        assert(ws);
        double weight=0.01;
        double nt(0);
        for (int i=0;i<this->dim_t;++i)
            nt += SQR(this->t[i]);
        double nr(0);
        for (int i=0;i<this->dim_r;++i)
            nr += SQR(this->r[i]);
        return nt + weight*nr;
    }

    Config normalize () const {
        assert(ws);

        Config nc = (*this);
        double n = norm();
        if (n > 0) {
            for (int i=0;i<nc.dim_t;++i)
                nc.t[i] = nc.t[i]/n;
            for (int i=0;i<nc.dim_r;++i)
                nc.r[i] = nc.r[i]/n;
        }
        return nc;
    }

    Config operator+ (const Config& rhs) const;
    Config operator- (const Config& rhs) const;
    Config operator* (const double rhs) const;
    Config operator/ (const double rhs) const;
    bool operator==(const Config& cfg) const;

    bool inbound () const;
    void round ();
    double eulerDiff (double e1, double e2) const;
    double translationRound (double t) const;
    double eulerRound (double e) const;

    friend std::ostream& operator<< (std::ostream& o, const Config& c) {
        o << "trans dim "  << c.dim_t << "\n";
        for(int i=0;i<c.dim_t;++i)
            o << c.t[i] << " ";
        o << "\n";
        o << "rotate dim " << c.dim_r << "\n";
        for(int i=0;i<c.dim_r;++i)
            o << c.r[i] << " ";
        o << "\n";
        return o;
    }

    static Config randomCfg (int dim_t, int dim_r);
    Config toPhysical (std::vector<double> &delta, std::vector<double> &ref) const;
    Config toParamtric (std::vector<double> &delta, std::vector<double> &ref) const;
    // distance in parametric space
    inline double distance (Config cfg) {
        //assert(ws);
        assert(dim_t == cfg.dim_t);
        double nt(0);
        for (int i=0;i<dim_t;++i)
            nt += SQR(t[i]-cfg.t[i]);
        return sqrt(nt);
    }


    int dim_t, dim_r;
    //all these values are between -1 and 1
    std::vector<double> t;
    std::vector<double> r;

    int source;
    bool ws; //true if the cfg is in parametric form. false if the cfg is in physical form.

private:
};

#endif // CONFIG_H
