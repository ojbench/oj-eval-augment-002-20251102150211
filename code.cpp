#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <string>

namespace sjtu {
class int2048 {
private:
  bool neg;                 // false for >=0, true for <0
  std::vector<int> d;       // least-significant block first
  static const int BASE = 10000;
  static const int WIDTH = 4;
  void trim();
  static int cmp_abs(const int2048 &a, const int2048 &b);
  static int2048 add_abs(const int2048 &a, const int2048 &b);
  static int2048 sub_abs(const int2048 &a, const int2048 &b); // |a| >= |b|
  static void divmod_abs(const int2048 &a, const int2048 &b, int2048 &q, int2048 &r);
public:
  int2048();
  int2048(long long);
  int2048(const std::string &);
  int2048(const int2048 &);
  void read(const std::string &);
  void print();
  int2048 &add(const int2048 &);
  friend int2048 add(int2048, const int2048 &);
  int2048 &minus(const int2048 &);
  friend int2048 minus(int2048, const int2048 &);
  int2048 operator+() const;
  int2048 operator-() const;
  int2048 &operator=(const int2048 &);
  int2048 &operator+=(const int2048 &);
  friend int2048 operator+(int2048, const int2048 &);
  int2048 &operator-=(const int2048 &);
  friend int2048 operator-(int2048, const int2048 &);
  int2048 &operator*=(const int2048 &);
  friend int2048 operator*(int2048, const int2048 &);
  int2048 &operator/=(const int2048 &);
  friend int2048 operator/(int2048, const int2048 &);
  int2048 &operator%=(const int2048 &);
  friend int2048 operator%(int2048, const int2048 &);
  friend std::istream &operator>>(std::istream &, int2048 &);
  friend std::ostream &operator<<(std::ostream &, const int2048 &);
  friend bool operator==(const int2048 &, const int2048 &);
  friend bool operator!=(const int2048 &, const int2048 &);
  friend bool operator<(const int2048 &, const int2048 &);
  friend bool operator>(const int2048 &, const int2048 &);
  friend bool operator<=(const int2048 &, const int2048 &);
  friend bool operator>=(const int2048 &, const int2048 &);
};
} // namespace sjtu

// Implementation
namespace sjtu {
using std::complex; using std::istream; using std::ostream; using std::string; using std::vector; using std::cout; using std::cin;
static const double PI=3.1415926535897932384626433832795;
static void fft(vector<complex<double>>& a,bool invert){int n=a.size(); static vector<int> rev; static vector<complex<double>> roots{0,1}; if((int)rev.size()!=n){int k=__builtin_ctz(n); rev.assign(n,0); for(int i=0;i<n;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(k-1));}
 if((int)roots.size()<n){int k=__builtin_ctz(roots.size()); roots.resize(n); while((1<<k)<n){double ang=2*PI/(1<<(k+1)); for(int i=1<<(k-1);i<(1<<k);++i){roots[2*i]=roots[i]; double a=ang*(2*i+1-(1<<k)); roots[2*i+1]=complex<double>(cos(a),sin(a));} ++k;}}
 for(int i=0;i<n;i++) if(i<rev[i]) std::swap(a[i],a[rev[i]]);
 for(int len=1;len<n;len<<=1){for(int i=0;i<n;i+=2*len){for(int j=0;j<len;j++){complex<double> u=a[i+j]; complex<double> v=a[i+j+len]*(invert?std::conj(roots[len+j]):roots[len+j]); a[i+j]=u+v; a[i+j+len]=u-v;}}} if(invert){for(int i=0;i<n;i++) a[i]/=n;}}
static vector<int> mul_conv_base(const vector<int>& A,const vector<int>& B){ if(A.empty()||B.empty()) return {}; int n1=A.size(), n2=B.size(); int n=1; while(n<n1+n2) n<<=1; vector<complex<double>> fa(n),fb(n); for(int i=0;i<n1;i++) fa[i]=A[i]; for(int i=0;i<n2;i++) fb[i]=B[i]; fft(fa,false); fft(fb,false); for(int i=0;i<n;i++) fa[i]*=fb[i]; fft(fa,true); vector<int> res(n); long long carry=0; for(int i=0;i<n;i++){ long long t=(long long)llround(fa[i].real())+carry; res[i]=int(t%10000); carry=t/10000;} while(carry){res.push_back(int(carry%10000)); carry/=10000;} while(res.size()>1&&res.back()==0) res.pop_back(); return res; }
void int2048::trim(){ while(!d.empty()&&d.back()==0) d.pop_back(); if(d.empty()) neg=false; }
int int2048::cmp_abs(const int2048& a,const int2048& b){ if(a.d.size()!=b.d.size()) return a.d.size()<b.d.size()?-1:1; for(int i=(int)a.d.size()-1;i>=0;--i) if(a.d[i]!=b.d[i]) return a.d[i]<b.d[i]?-1:1; return 0; }
int2048 int2048::add_abs(const int2048& a,const int2048& b){ int2048 r; r.neg=false; int n=std::max(a.d.size(),b.d.size()); r.d.assign(n+1,0); int carry=0; for(int i=0;i<n;i++){ int x=i<(int)a.d.size()?a.d[i]:0; int y=i<(int)b.d.size()?b.d[i]:0; int s=x+y+carry; r.d[i]=s%BASE; carry=s/BASE;} r.d[n]=carry; r.trim(); return r; }
int2048 int2048::sub_abs(const int2048& a,const int2048& b){ int2048 r; r.neg=false; int n=a.d.size(); r.d.assign(n,0); int carry=0; for(int i=0;i<n;i++){ int x=a.d[i]-carry-(i<(int)b.d.size()?b.d[i]:0); if(x<0){x+=BASE; carry=1;} else carry=0; r.d[i]=x;} r.trim(); return r; }
int2048::int2048(){neg=false;}
int2048::int2048(long long x){neg=false; d.clear(); if(x<0){neg=true; x=-x;} if(x==0){} else{ while(x){d.push_back((int)(x%BASE)); x/=BASE;}}}
int2048::int2048(const string& s){read(s);} 
int2048::int2048(const int2048& o){neg=o.neg; d=o.d;}
void int2048::read(const string& s){neg=false; d.clear(); int i=0; while(i<(int)s.size()&&isspace((unsigned char)s[i]))++i; bool negs=false; if(i<(int)s.size()&&(s[i]=='+'||s[i]=='-')){negs=s[i]=='-'; ++i;} while(i<(int)s.size()&&s[i]=='0') ++i; vector<int> tmp; for(int j=s.size()-1;j>=i;j-=WIDTH){ int x=0; int l=std::max(i,j-WIDTH+1); for(int k=l;k<=j;k++) x=x*10+(s[k]-'0'); tmp.push_back(x);} d=tmp; trim(); if(!d.empty()) neg=negs; }
void int2048::print(){ if(d.empty()){cout<<0; return;} if(neg) cout<<'-'; cout<<d.back(); for(int i=(int)d.size()-2;i>=0;--i){ int x=d[i]; if(x<1000) cout<<'0'; if(x<100) cout<<'0'; if(x<10) cout<<'0'; cout<<x;} }
int2048& int2048::add(const int2048& b){ if(b.d.empty()) return *this; if(d.empty()){*this=b; return *this;} bool A=neg, B=b.neg; if(A==B){ int2048 r=add_abs(*this,b); r.neg=A; *this=r; } else { int c=cmp_abs(*this,b); if(c==0){ d.clear(); neg=false; } else if(c>0){ int2048 r=sub_abs(*this,b); r.neg=A; *this=r; } else { int2048 r=sub_abs(b,*this); r.neg=B; *this=r; } } return *this; }
int2048 add(int2048 a,const int2048& b){ return a.add(b); }
int2048& int2048::minus(const int2048& b){ int2048 nb=b; nb.neg=!b.neg; return add(nb); }
int2048 minus(int2048 a,const int2048& b){ return a.minus(b); }
int2048 int2048::operator+() const{ return *this; }
int2048 int2048::operator-() const{ int2048 r=*this; if(!r.d.empty()) r.neg=!r.neg; return r; }
int2048& int2048::operator=(const int2048& o){ if(this!=&o){neg=o.neg; d=o.d;} return *this; }
int2048& int2048::operator+=(const int2048& b){ return add(b);} 
int2048 operator+(int2048 a,const int2048& b){ return a+=b; }
int2048& int2048::operator-=(const int2048& b){ return minus(b);} 
int2048 operator-(int2048 a,const int2048& b){ return a-=b; }
int2048& int2048::operator*=(const int2048& b){ int n=d.size(), m=b.d.size(); if(!n||!m){ d.clear(); neg=false; return *this;} int2048 r; r.neg=false; if(1ll*n*m<=4096){ r.d.assign(n+m,0); for(int i=0;i<n;i++){ long long carry=0; for(int j=0;j<m;j++){ long long cur=r.d[i+j]+1ll*d[i]*b.d[j]+carry; r.d[i+j]=int(cur%BASE); carry=cur/BASE;} int k=i+m; while(carry){ long long cur=r.d[k]+carry; r.d[k]=int(cur%BASE); carry=cur/BASE; ++k;} } } else { vector<int> res=mul_conv_base(d,b.d); r.d.swap(res);} r.trim(); r.neg = (!r.d.empty()) && (neg^b.neg); *this=r; return *this; }
int2048 operator*(int2048 a,const int2048& b){ return a*=b; }
void int2048::divmod_abs(const int2048& a,const int2048& b,int2048& q,int2048& r){ q.d.assign(a.d.size(),0); q.neg=false; r.neg=false; r.d.clear(); if(b.d.empty()){ q.d.clear(); r.d.clear(); return;} int n=a.d.size(); for(int i=n-1;i>=0;--i){ if(!r.d.empty()||a.d[i]) r.d.insert(r.d.begin(),a.d[i]); r.trim(); int lo=0,hi=BASE-1,ans=0; auto mul_small=[&](const int2048& x,int m){ int2048 t; t.neg=false; if(m==0||x.d.empty()) return t; t.d.assign(x.d.size()+1,0); long long carry=0; for(size_t k=0;k<x.d.size();++k){ long long cur=carry+1ll*x.d[k]*m; t.d[k]=int(cur%BASE); carry=cur/BASE;} if(carry) t.d[x.d.size()]=int(carry); t.trim(); return t; }; while(lo<=hi){ int mid=(lo+hi)>>1; int2048 t=mul_small(b,mid); int c=cmp_abs(t,r); if(c<=0){ ans=mid; lo=mid+1;} else hi=mid-1; } if(ans){ int2048 t=mul_small(b,ans); r=sub_abs(r,t); } q.d[i]=ans; } q.trim(); r.trim(); }
int2048& int2048::operator/=(const int2048& b){ if(b.d.empty()){ *this=int2048(); return *this;} if(d.empty()) return *this; int2048 a=*this, bb=b; a.neg=false; bb.neg=false; int2048 q,r; divmod_abs(a,bb,q,r); bool negres=(this->neg^b.neg); if(negres && !r.d.empty()){ q=add_abs(q,int2048(1)); r=sub_abs(bb,r); } q.neg=negres && !q.d.empty(); *this=q; return *this; }
int2048 operator/(int2048 a,const int2048& b){ return a/=b; }
int2048& int2048::operator%=(const int2048& b){ if(b.d.empty()){ d.clear(); neg=false; return *this;} if(d.empty()) return *this; int2048 a=*this, bb=b; a.neg=false; bb.neg=false; int2048 q,r; divmod_abs(a,bb,q,r); if((this->neg^b.neg) && !r.d.empty()){ r=sub_abs(bb,r); } r.neg=b.neg?(!r.d.empty()):false; *this=r; return *this; }
int2048 operator%(int2048 a,const int2048& b){ return a%=b; }
istream& operator>>(istream& is,int2048& x){ string s; is>>s; x.read(s); return is; }
ostream& operator<<(ostream& os,const int2048& x){ if(x.d.empty()) return os<<0; if(x.neg) os<<'-'; os<<x.d.back(); for(int i=(int)x.d.size()-2;i>=0;--i){ int v=x.d[i]; if(v<1000) os<<'0'; if(v<100) os<<'0'; if(v<10) os<<'0'; os<<v;} return os; }
bool operator==(const int2048& a,const int2048& b){ return a.neg==b.neg && a.d==b.d; }
bool operator!=(const int2048& a,const int2048& b){ return !(a==b); }
bool operator<(const int2048& a,const int2048& b){ if(a.neg!=b.neg) return a.neg; int c=int2048::cmp_abs(a,b); return a.neg?c>0:c<0; }
bool operator>(const int2048& a,const int2048& b){ return b<a; }
bool operator<=(const int2048& a,const int2048& b){ return !(b<a); }
bool operator>=(const int2048& a,const int2048& b){ return !(a<b); }
} // namespace sjtu

