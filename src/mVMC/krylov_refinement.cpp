/*
 * Exact binary64-input algebra and directed-rounding exponential enclosures
 * for unresolved GC final-state decisions. GPL-3.0-or-later, as mVMC.
 */
#include <gmpxx.h>
#include <mpfr.h>
#include <algorithm>
#include <cmath>
#include <complex.h>
#include <cstring>
#include <limits>
#include <map>
#include <new>
#include <string>
#include <utility>
#include <vector>

/* GCC/Clang C complex ABI, shared with the C sampler. */
#define complex __complex__
extern "C" {
#include "krylov_refinement.h"
}
#undef complex

namespace {
typedef __complex__ double CComplex;
struct Failure { MVMCKrylovStatus status; };
void require(bool condition, MVMCKrylovStatus status =
             MVMC_KRYLOV_STATUS_INVALID_ARGUMENT) {
  if (!condition) throw Failure{status};
}
bool finite(CComplex z) {
  return std::isfinite(__real__ z) && std::isfinite(__imag__ z);
}
struct QComplex {
  mpq_class r, i;
  QComplex() {}
  explicit QComplex(long x) : r(x) {}
  explicit QComplex(CComplex z) {
    require(finite(z), MVMC_KRYLOV_STATUS_NONFINITE);
    mpq_set_d(r.get_mpq_t(), __real__ z);
    mpq_set_d(i.get_mpq_t(), __imag__ z);
  }
  QComplex(const mpq_class &a, const mpq_class &b) : r(a), i(b) {}
  bool zero() const { return r == 0 && i == 0; }
};
QComplex operator+(const QComplex &a, const QComplex &b) {
  return QComplex(a.r+b.r, a.i+b.i);
}
QComplex operator-(const QComplex &a, const QComplex &b) {
  return QComplex(a.r-b.r, a.i-b.i);
}
QComplex operator-(const QComplex &a) { return QComplex(-a.r, -a.i); }
QComplex operator*(const QComplex &a, const QComplex &b) {
  return QComplex(a.r*b.r-a.i*b.i, a.r*b.i+a.i*b.r);
}
QComplex operator/(const QComplex &a, const QComplex &b) {
  mpq_class d = b.r*b.r+b.i*b.i;
  require(d != 0, MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE);
  return QComplex((a.r*b.r+a.i*b.i)/d, (a.i*b.r-a.r*b.i)/d);
}
void check(const mpq_class &x, size_t bits) {
  require(mpz_sizeinbase(mpq_numref(x.get_mpq_t()), 2) <= bits &&
          mpz_sizeinbase(mpq_denref(x.get_mpq_t()), 2) <= bits,
          MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
}
void check(const QComplex &x, size_t bits) {
  check(x.r, bits); check(x.i, bits);
}
typedef std::map<mpq_class, QComplex> Expression;
typedef std::vector<uint64_t> Words;
const MVMCKrylovRefinementLimits defaults = {
  16384, 8192, 8192, 65536, UINT64_C(2000000), 4096
};

/* Every interval operation includes conversion and arithmetic rounding.
 * No tolerance is used to declare a mathematical zero. */
struct Real {
  mpfr_t x;
  explicit Real(unsigned p) { mpfr_init2(x, p); mpfr_set_zero(x, 1); }
  Real(const Real &a) {
    mpfr_init2(x, mpfr_get_prec(a.x)); mpfr_set(x, a.x, MPFR_RNDN);
  }
  Real &operator=(const Real &a) {
    mpfr_set_prec(x, mpfr_get_prec(a.x)); mpfr_set(x, a.x, MPFR_RNDN);
    return *this;
  }
  ~Real() { mpfr_clear(x); }
};
struct Interval {
  Real lo, hi;
  explicit Interval(unsigned p) : lo(p), hi(p) {}
  Interval(const mpq_class &q, unsigned p) : lo(p), hi(p) {
    mpfr_set_q(lo.x, q.get_mpq_t(), MPFR_RNDD);
    mpfr_set_q(hi.x, q.get_mpq_t(), MPFR_RNDU);
  }
  unsigned precision() const { return (unsigned)mpfr_get_prec(lo.x); }
};
Interval add(const Interval &a, const Interval &b) {
  Interval c(a.precision());
  mpfr_add(c.lo.x,a.lo.x,b.lo.x,MPFR_RNDD);
  mpfr_add(c.hi.x,a.hi.x,b.hi.x,MPFR_RNDU); return c;
}
Interval neg(const Interval &a) {
  Interval c(a.precision());
  mpfr_neg(c.lo.x,a.hi.x,MPFR_RNDN);
  mpfr_neg(c.hi.x,a.lo.x,MPFR_RNDN); return c;
}
Interval mul(const Interval &a, const Interval &b) {
  Interval c(a.precision()); Real l(a.precision()), h(a.precision());
  mpfr_set_inf(c.lo.x,1); mpfr_set_inf(c.hi.x,-1);
  for (int j=0;j<2;++j) for (int k=0;k<2;++k) {
    mpfr_srcptr x=j ? a.hi.x : a.lo.x, y=k ? b.hi.x : b.lo.x;
    mpfr_mul(l.x,x,y,MPFR_RNDD); mpfr_mul(h.x,x,y,MPFR_RNDU);
    require(!mpfr_nan_p(l.x) && !mpfr_nan_p(h.x),
            MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
    if(mpfr_less_p(l.x,c.lo.x)) mpfr_set(c.lo.x,l.x,MPFR_RNDN);
    if(mpfr_greater_p(h.x,c.hi.x)) mpfr_set(c.hi.x,h.x,MPFR_RNDN);
  }
  return c;
}
Interval square(const Interval &a) {
  Interval c=mul(a,a);
  if(mpfr_sgn(a.lo.x)<=0 && mpfr_sgn(a.hi.x)>=0)
    mpfr_set_zero(c.lo.x,1);
  return c;
}
Interval divide_positive(const Interval &a, const Interval &b) {
  require(mpfr_sgn(b.lo.x)>0, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
  Interval inv(b.precision());
  mpfr_ui_div(inv.lo.x,1,b.hi.x,MPFR_RNDD);
  mpfr_ui_div(inv.hi.x,1,b.lo.x,MPFR_RNDU); return mul(a,inv);
}
struct ComplexInterval {
  Interval r,i;
  explicit ComplexInterval(unsigned p): r(p),i(p) {}
};
Interval norm(const ComplexInterval &z) {
  return add(square(z.r),square(z.i));
}
ComplexInterval evaluate(const Expression &e, const mpq_class &shift,
                         unsigned p) {
  ComplexInterval z(p);
  for (const auto &t : e) {
    Interval g(t.first-shift,p), exponential(p);
    mpfr_exp(exponential.lo.x,g.lo.x,MPFR_RNDD);
    mpfr_exp(exponential.hi.x,g.hi.x,MPFR_RNDU);
    z.r=add(z.r,mul(Interval(t.second.r,p),exponential));
    z.i=add(z.i,mul(Interval(t.second.i,p),exponential));
  }
  return z;
}
ComplexInterval quotient(const ComplexInterval &a, const ComplexInterval &b) {
  ComplexInterval c(a.r.precision()); Interval d=norm(b);
  c.r=divide_positive(add(mul(a.r,b.r),mul(a.i,b.i)),d);
  c.i=divide_positive(add(mul(a.i,b.r),neg(mul(a.r,b.i))),d);
  return c;
}
bool narrow(const ComplexInterval &z) {
  const unsigned p=z.r.precision();
  Interval n=norm(z); Real width(p), w(p), scale(p);
  if (!mpfr_number_p(n.hi.x) || mpfr_sgn(n.lo.x)<=0) return false;
  mpfr_sub(width.x,z.r.hi.x,z.r.lo.x,MPFR_RNDU);
  mpfr_sub(w.x,z.i.hi.x,z.i.lo.x,MPFR_RNDU);
  mpfr_add(width.x,width.x,w.x,MPFR_RNDU);
  mpfr_sqrt(scale.x,n.lo.x,MPFR_RNDD);
  mpfr_mul_2si(scale.x,scale.x,-80,MPFR_RNDD);
  return mpfr_lessequal_p(width.x,scale.x);
}
MVMCScaledComplex export_value(const ComplexInterval &z,
                               const mpq_class &shift, bool exact_zero) {
  MVMCScaledComplex out;
  if (exact_zero) {
    mvmc_scaled_complex_make_exact_zero(&out); return out;
  }
  require(narrow(z), MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
  unsigned p=z.r.precision();
  Real r(p),i(p),a(p),b(p),logabs(p),g(p);
  mpfr_add(r.x,z.r.lo.x,z.r.hi.x,MPFR_RNDN);
  mpfr_div_2ui(r.x,r.x,1,MPFR_RNDN);
  mpfr_add(i.x,z.i.lo.x,z.i.hi.x,MPFR_RNDN);
  mpfr_div_2ui(i.x,i.x,1,MPFR_RNDN);
  mpfr_hypot(a.x,r.x,i.x,MPFR_RNDN);
  mpfr_div(b.x,r.x,a.x,MPFR_RNDN);
  CComplex phase; __real__ phase=mpfr_get_d(b.x,MPFR_RNDN);
  mpfr_div(b.x,i.x,a.x,MPFR_RNDN);
  __imag__ phase=mpfr_get_d(b.x,MPFR_RNDN);
  mpfr_log(logabs.x,a.x,MPFR_RNDN);
  mpfr_set_q(g.x,shift.get_mpq_t(),MPFR_RNDN);
  mpfr_add(logabs.x,logabs.x,g.x,MPFR_RNDN);
  double logvalue=mpfr_get_d(logabs.x,MPFR_RNDN);
  require(std::isfinite(logvalue),MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
  /* Export is the ordinary phase/log representation, not the certificate.
   * Include log-coordinate rounding in its conservative error estimate. */
  double relative=64*std::numeric_limits<double>::epsilon()*
                  (1+std::fabs(logvalue));
  require(relative < 0.01,MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
  require(mvmc_scaled_complex_make_finite(phase,logvalue,
            logvalue+std::log(relative),&out)==MVMC_PFAFFIAN_STATUS_OK &&
          out.state==MVMC_SCALED_COMPLEX_FINITE_NONZERO,
          MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
  return out;
}
} // namespace

struct MVMCRefinedBase {
  QComplex coefficient;
  mpq_class exponent;
  const MVMCKrylovRefinementLimits *limits;
};

extern "C" MVMCKrylovStatus mvmc_refined_base_build(
    MVMCRefinedBase *result,const CComplex *slater,size_t norb,
    const int *occupied,size_t n,const CComplex *weights,size_t nqp,
    int support_zero,const CComplex *parameters,const int64_t *counts,
    size_t nproj) {
  try {
    require(result && result->limits && slater && weights && nqp &&
            norb && n<=norb && n%2==0 && (!n || occupied) &&
            (!nproj || (parameters && counts)));
    const auto &limits=*result->limits;
    require(norb<=SIZE_MAX/norb && nqp<=SIZE_MAX/(norb*norb) &&
            n<=limits.max_matrix_entries/(n ? n : 1),
            MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
    for (size_t k=0;k<n;++k)
      require(occupied[k]>=0 && (size_t)occupied[k]<norb &&
              (!k || occupied[k]>occupied[k-1]));
    /* Validate the full immutable input, including support-zero cases. */
    for (size_t k=0;k<nqp*norb*norb;++k)
      require(finite(slater[k]),MVMC_KRYLOV_STATUS_NONFINITE);
    for (size_t k=0;k<nqp;++k)
      require(finite(weights[k]),MVMC_KRYLOV_STATUS_NONFINITE);
    mpq_class g;
    for (size_t k=0;k<nproj;++k) {
      require(finite(parameters[k]),MVMC_KRYLOV_STATUS_NONFINITE);
      mpq_class parameter; mpq_set_d(parameter.get_mpq_t(),__real__ parameters[k]);
      /* Portable exact signed 64-bit conversion (also on LLP64). */
      mpz_class count(std::to_string(counts[k]));
      g += parameter*count; check(g,limits.max_rational_bits);
    }
    QComplex total;
    if (!support_zero) for (size_t qp=0;qp<nqp;++qp) {
      std::vector<QComplex> a(n*n);
      for (size_t j=0;j<n;++j) for (size_t k=0;k<n;++k)
        a[j*n+k]=QComplex(slater[qp*norb*norb+
                                (size_t)occupied[j]*norb+occupied[k]]);
      for (size_t j=0;j<n;++j) {
        require(a[j*n+j].zero());
        for (size_t k=j+1;k<n;++k)
          require((a[j*n+k]+a[k*n+j]).zero());
      }
      QComplex pf(1L);
      for (size_t j=0;j<n;j+=2) {
        size_t pivot=j+1;
        while(pivot<n && a[j*n+pivot].zero()) ++pivot;
        if(pivot==n) { pf=QComplex(); break; }
        if(pivot!=j+1) {
          for(size_t k=0;k<n;++k) std::swap(a[(j+1)*n+k],a[pivot*n+k]);
          for(size_t k=0;k<n;++k) std::swap(a[k*n+j+1],a[k*n+pivot]);
          pf=-pf;
        }
        QComplex d=a[j*n+j+1]; pf=pf*d;
        check(pf,limits.max_rational_bits);
        for(size_t k=j+2;k<n;++k) for(size_t l=k+1;l<n;++l) {
          a[k*n+l]=a[k*n+l]-
            (a[j*n+k]*a[(j+1)*n+l]-a[j*n+l]*a[(j+1)*n+k])/d;
          check(a[k*n+l],limits.max_rational_bits);
          a[l*n+k]=-a[k*n+l];
        }
      }
      total=total+QComplex(weights[qp])*pf;
      check(total,limits.max_rational_bits);
    }
    result->coefficient=total; result->exponent=g;
    return MVMC_KRYLOV_STATUS_OK;
  } catch(const Failure &f) { return f.status; }
    catch(const std::bad_alloc &) { return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE; }
}

namespace {
struct Engine {
  const MVMCKrylovFockModel &model;
  MVMCRefinedBaseCallback amplitude;
  void *context;
  const MVMCKrylovRefinementLimits &limits;
  uint64_t calls=0;
  size_t stored_terms=0;
  std::map<std::pair<int,Words>,Expression> cache;
  Engine(const MVMCKrylovFockModel &m, MVMCRefinedBaseCallback a, void *c,
         const MVMCKrylovRefinementLimits &l) :
    model(m),amplitude(a),context(c),limits(l) {}
  void accumulate(Expression &a,const Expression &b,const QComplex &factor) {
    if(factor.zero()) return;
    for(const auto &t:b) {
      QComplex value=a[t.first]+factor*t.second;
      check(value,limits.max_rational_bits);
      if(value.zero()) a.erase(t.first); else a[t.first]=value;
      require(a.size()<=limits.max_expression_terms,
              MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
    }
  }
  Expression evaluate(const Words &words,int depth) {
    require(calls<limits.max_node_calls,MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
    ++calls;
    auto key=std::make_pair(depth,words);
    auto found=cache.find(key); if(found!=cache.end()) return found->second;
    Expression value;
    if(depth==0) {
      MVMCRefinedBase base; base.limits=&limits;
      auto status=amplitude(words.data(),words.size(),context,&base);
      require(status==MVMC_KRYLOV_STATUS_OK,status);
      if(!base.coefficient.zero()) value[base.exponent]=base.coefficient;
    } else {
      for(size_t t=0;t<model.term_count;++t) {
        const auto &term=model.terms[t]; Words next=words;
        int sign=1; bool applied=true;
        for(size_t k=term.operator_count;k>0;--k) {
          const auto &op=model.operators[term.operator_offset+k-1];
          bool occupied=(next[op.orbital/64]>>(op.orbital%64))&1;
          if(occupied==(op.kind==MVMC_KRYLOV_FERMION_CREATE)) {
            applied=false; break;
          }
          unsigned parity=0;
          for(size_t l=0;l<op.orbital/64;++l)
            parity^=(unsigned)__builtin_popcountll(next[l])&1;
          uint64_t mask=(UINT64_C(1)<<(op.orbital%64))-1;
          parity^=(unsigned)__builtin_popcountll(next[op.orbital/64]&mask)&1;
          if(parity) sign=-sign;
          next[op.orbital/64]^=UINT64_C(1)<<(op.orbital%64);
        }
        if(applied) {
          /* Hermiticity turns the outgoing matrix element into this row. */
          QComplex coefficient(term.coefficient);
          coefficient.i=-coefficient.i;
          if(sign<0) coefficient=-coefficient;
          if(!coefficient.zero())
            accumulate(value,evaluate(next,depth-1),coefficient);
        }
      }
    }
    require(cache.size()<limits.max_cache_entries &&
            value.size()<=limits.max_expression_terms-stored_terms,
            MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
    stored_terms+=value.size(); cache.emplace(key,value); return value;
  }
};
mpq_class maximum_exponent(const Expression &a,const Expression &b) {
  require(!a.empty(),MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE);
  mpq_class g=a.rbegin()->first;
  if(!b.empty() && b.rbegin()->first>g) g=b.rbegin()->first;
  return g;
}
} // namespace

extern "C" MVMCKrylovStatus mvmc_krylov_refine_final(
    const MVMCKrylovFockModel *model,const MVMCKrylovFinalStatePolicy *policy,
    const uint64_t *current,const uint64_t *proposal,size_t nw,double uniform,
    MVMCRefinedBaseCallback amplitude,void *context,
    const MVMCKrylovRefinementLimits *requested,
    MVMCKrylovRefinementResult *result) {
  try {
    require(model && policy && current && proposal && nw && amplitude && result &&
            std::isfinite(uniform) && uniform>=0 && uniform<1 &&
            policy->order>=1 && policy->order<=2 &&
            policy->policy_hash &&
            policy->policy_hash==mvmc_krylov_final_state_policy_hash(policy));
    require(model->hermitian, MVMC_KRYLOV_STATUS_NON_HERMITIAN_MODEL);
    require(model->site_count && model->site_count<=SIZE_MAX/2 &&
            nw==(2*model->site_count/64+(2*model->site_count%64!=0)) &&
            (!model->term_count || model->terms) &&
            (!model->operator_count || model->operators));
    require(model->domain==MVMC_KRYLOV_DOMAIN_GC_BALANCED ||
            model->domain==MVMC_KRYLOV_DOMAIN_GC_EVEN,
            MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL);
    MVMCKrylovStatus validation=mvmc_krylov_fock_validate(model,current,nw);
    require(validation==MVMC_KRYLOV_STATUS_OK,validation);
    validation=mvmc_krylov_fock_validate(model,proposal,nw);
    require(validation==MVMC_KRYLOV_STATUS_OK,validation);
    for(size_t k=0;k<model->term_count;++k) {
      const auto &t=model->terms[k];
      require(finite(t.coefficient),MVMC_KRYLOV_STATUS_NONFINITE);
      require(t.operator_offset<=model->operator_count &&
              t.operator_count<=model->operator_count-t.operator_offset);
    }
    for(size_t k=0;k<model->operator_count;++k)
      require(model->operators[k].orbital<2*model->site_count &&
              (model->operators[k].kind==MVMC_KRYLOV_FERMION_CREATE ||
               model->operators[k].kind==MVMC_KRYLOV_FERMION_ANNIHILATE));
    const auto &limits=requested ? *requested : defaults;
    require(limits.max_matrix_entries && limits.max_cache_entries &&
            limits.max_expression_terms && limits.max_rational_bits &&
            limits.max_node_calls && limits.max_precision>=128 &&
            limits.max_precision<=65536);
    Engine engine(*model,amplitude,context,limits);
    Words x(current,current+nw), y(proposal,proposal+nw);
    Expression phi_x,phi_y,hphi_y;
    std::vector<Expression> basis;
    for(int k=0;k<=policy->order;++k) basis.push_back(engine.evaluate(y,k));
    for(int k=0;k<=policy->order;++k) {
      QComplex c(policy->coefficient[k]);
      engine.accumulate(phi_y,basis[k],c);
    }
    MVMCKrylovRefinementResult candidate={};
    candidate.node_calls=engine.calls;
    if(phi_y.empty()) {
      candidate.valid=1; candidate.exact_zero_proposal=1;
      *result=candidate; return MVMC_KRYLOV_STATUS_OK;
    }
    for(int k=0;k<=policy->order;++k)
      engine.accumulate(phi_x,engine.evaluate(x,k),QComplex(policy->coefficient[k]));
    require(!phi_x.empty(),MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE);
    candidate.node_calls=engine.calls;
    const mpq_class shift=maximum_exponent(phi_x,phi_y);
    mpq_class u; mpq_set_d(u.get_mpq_t(),uniform);
    bool rational_accept=false;
    if(phi_x.size()==1 && phi_y.size()==1 &&
       phi_x.begin()->first==phi_y.begin()->first) {
      const auto &xq=phi_x.begin()->second, &yq=phi_y.begin()->second;
      mpq_class nx=xq.r*xq.r+xq.i*xq.i, ny=yq.r*yq.r+yq.i*yq.i;
      if(ny<=u*nx) {
        candidate.valid=1; *result=candidate;
        return MVMC_KRYLOV_STATUS_OK;
      }
      rational_accept=true;
    }
    for(unsigned p=128;p<=limits.max_precision;p*=2) {
      ComplexInterval a=evaluate(phi_x,shift,p), b=evaluate(phi_y,shift,p);
      Interval nx=norm(a), ny=norm(b), threshold=mul(Interval(u,p),nx);
      if(mpfr_sgn(nx.lo.x)<=0) continue;
      if(mpfr_lessequal_p(ny.hi.x,threshold.lo.x)) {
        candidate.valid=1; candidate.precision=p;
        *result=candidate; return MVMC_KRYLOV_STATUS_OK;
      }
      if(!rational_accept && !mpfr_greater_p(ny.lo.x,threshold.hi.x)) continue;
      /* A decision can be certified before energy is exportable. Refine
       * further using the same draw; never silently accept bad energy. */
      if(basis.size()==(size_t)policy->order+1) {
        basis.push_back(engine.evaluate(y,policy->order+1));
        for(int k=0;k<=policy->order;++k)
          engine.accumulate(hphi_y,basis[k+1],QComplex(policy->coefficient[k]));
        candidate.node_calls=engine.calls;
      }
      const mpq_class sy=maximum_exponent(phi_y,hphi_y);
      b=evaluate(phi_y,sy,p);
      ComplexInterval h=evaluate(hphi_y,sy,p);
      if(!narrow(b) || (!hphi_y.empty() && !narrow(h))) continue;
      ComplexInterval energy=quotient(h,b);
      if(!hphi_y.empty() && !narrow(energy)) continue;
      bool ready=true;
      for(const auto &v:basis)
        if(!v.empty() && !narrow(evaluate(v,v.rbegin()->first,p))) ready=false;
      if(!ready) continue;
      candidate.valid=1; candidate.accepted=1; candidate.precision=p;
      auto &krylov=candidate.krylov; auto &e=candidate.evaluation;
      krylov.valid=1; krylov.status=MVMC_KRYLOV_STATUS_OK;
      krylov.evaluated_order=policy->order+1;
      for(size_t k=0;k<basis.size();++k) {
        const auto &v=basis[k]; mpq_class s=v.empty()?mpq_class(0):v.rbegin()->first;
        krylov.value[k]=export_value(evaluate(v,s,p),s,v.empty());
        if(v.empty()) ++e.basis_exact_zero_count; else ++e.basis_finite_count;
      }
      e.valid=1; e.status=MVMC_KRYLOV_STATUS_OK; e.order=policy->order;
      e.policy_hash=policy->policy_hash; e.sampleable=1; e.local_energy_available=1;
      for(int k=0;k<=policy->order;++k)
        if(QComplex(policy->coefficient[k]).zero()) ++e.coefficient_zero_count;
      e.amplitude=export_value(b,sy,false);
      e.hamiltonian_amplitude=export_value(h,sy,hphi_y.empty());
      e.local_energy=export_value(energy,mpq_class(0),hphi_y.empty());
      e.log_weight=2*e.amplitude.log_abs;
      require(std::isfinite(e.log_weight),MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
      *result=candidate; return MVMC_KRYLOV_STATUS_OK;
    }
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  } catch(const Failure &f) { return f.status; }
    catch(const std::bad_alloc &) { return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE; }
}
