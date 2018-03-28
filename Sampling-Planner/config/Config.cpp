/*
 *  Config.cpp
 *
 *  Created on: Mar. 22, 2018
 *      Author: Ching-Hsiang Hsu
 *
 */

#include "Config.hpp"

Config::Config () {
    dim_t=dim_r=0;
    t.clear();
    r.clear();
    source=-1;
    ws=true;
}

Config::Config (const Config &c) {
    copy(c.dim_t, c.t, this->dim_t, this->t);
    copy(c.dim_r, c.r, this->dim_r, this->r);
    source=c.source;
    ws=c.ws;
}

Config::Config (int dim_t, int dim_r, int source, bool ws) {
    init(dim_t, t, this->dim_t, this->t);
    init(dim_r, r, this->dim_r, this->r);

    this->source = source;
    this->ws = ws;
}

Config::Config (int dim_t, std::vector<double> t, int dim_r, std::vector<double> r, int source, bool ws) {
    copy(dim_t, t, this->dim_t, this->t);
    copy(dim_r, r, this->dim_r, this->r);

    this->source = source;
    this->ws = ws;
}

Config::~Config () {}

Config Config::operator+ (const Config& rhs) const {
    assert(ws);
    assert(rhs.ws);
    assert(this->dim_t == rhs.dim_t && this->dim_r == rhs.dim_r);

    Config r = (*this);
    for (int i=0;i<r.dim_t;++i)
        r.t[i] = r.t[i]+rhs.t[i];
    for (int i=0;i<r.dim_r;++i)
        r.r[i] = r.r[i]+rhs.r[i];
    r.round();
    return r;
}

Config Config::operator- (const Config& rhs) const {
    assert(ws);
    assert(rhs.ws);
    assert(dim_t == rhs.dim_t && dim_r == rhs.dim_r);

    Config r = (*this);
    for (int i=0;i<r.dim_t;++i)
        r.t[i] = r.t[i]-rhs.t[i];
    for (int i=0;i<r.dim_r;++i)
        r.r[i] = eulerDiff(r.r[i], rhs.r[i]);
    r.round();
    return r;
}

Config Config::operator* (const double rhs) const {
    assert(ws);

    Config r = (*this);
    for (int i=0;i<r.dim_t;++i)
        r.t[i] = r.t[i]*rhs;
    for (int i=0;i<r.dim_r;++i)
        r.r[i] = r.r[i]*rhs;
    r.round();
    return r;
}

Config Config::operator/ (const double rhs) const {
    Config r = (*this);
    if(rhs != 0) {
        double s = 1.0 / rhs;
        r = r*s;
    }
    return r;
}

//same config
bool Config::operator== (const Config& cfg) const {
    if (this->dim_t != cfg.dim_t) return false;
    if (this->dim_r != cfg.dim_r) return false;
    for(int i=0;i<this->dim_t;++i)
        if(this->t[i] != cfg.t[i]) return false;
    for(int i=0;i<this->dim_r;++i)
        if(this->r[i] != cfg.r[i]) return false;
    return true;
}

bool Config::inbound () const {
    assert(ws);

    for (int i=0;i<this->dim_t;++i)
        if(this->t[i]>1 || this->t[i]<-1) return false;
    for (int i=0;i<this->dim_r;++i)
        if(this->r[i]>1 || this->r[i]<-1) return false;
    return true;
}

void Config::round () {
    assert(ws);
    //for(int i=0;i<this->dim_r;++i)
    //this->t[i]=translationRound(this->t[i]);
    for(int i=0;i<dim_r;++i)
        r[i]=eulerRound(r[i]);
}

//compute e1 - e2
//e1 and e2 should be between -1 and 1
double Config::eulerDiff (double e1, double e2) const {
    double diff1 = (e1-e2);
    double diff2 = (e1<e2)?(e1+2-e2):(e1-e2-2);
    double d = (fabs(diff1)<fabs(diff2))?diff1:diff2;
    assert (d<=1&&d>=-1); //the result must be between -1 and 1
    return d;
}

//around e between -1 and 1
double Config::translationRound (double t) const {
    if (t>1) return 1;
    if (t<-1) return -1;
    return t;
}

//around e between -1 and 1
double Config::eulerRound (double e) const {
    return fmod((e+1),2)-1;
}

Config Config::randomCfg (int dim_t, int dim_r) {
    Config c(dim_t, dim_r);
    for (int i=0;i<dim_t;++i)
        c.t[i]=drand48()+drand48()-1;
    for (int i=0;i<dim_r;++i)
        c.r[i]=drand48()+drand48()-1;
    return c;
}

Config Config::toPhysical (std::vector<double> &delta, std::vector<double> &ref) const {
    assert(ws);

    Config nc = (*this);
    for (int i=0;i<nc.dim_t;++i)
        nc.t[i] = nc.t[i]*delta[i]*0.5+ref[i];
    for (int i=0;i<nc.dim_r;++i)
        nc.r[i] = nc.r[i]*180.0;
    nc.ws=false;
    nc.source=source;
    return nc;
}

Config Config::toParamtric (std::vector<double> &delta, std::vector<double> &ref) const {
    assert(!ws);

    Config nc = (*this);
    for (int i=0;i<nc.dim_t;++i)
        nc.t[i] = (nc.t[i]-ref[i])/(delta[i]*0.5);
    for (int i=0;i<nc.dim_r;++i)
        nc.r[i]=nc.eulerRound(nc.r[i]/180.0);
    nc.ws=true;
    nc.source=source;
    return nc;
}

