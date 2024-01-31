
#include <iostream>
#include <cstdlib>

#include <algorithm>
#include <vector>
#include <map>
#include <list>
#include <set>
#include <vector>
namespace linearExpr {
   using namespace std;
   const char names[27] = "xyzuvwabcdefghijklmnopqrst";

   template<class T>
   struct Eps {
      static constexpr T eps = 0;
   };
   template<>
   struct Eps<double> {
      static constexpr double eps = 1e-15;
   };
   template<>
   struct Eps<float> {
      static constexpr float eps = 1e-6f;
   };

   template<class T>
   bool ceq(T a, T b) {
      if (abs(a - b) < abs(a * Eps<T>::eps))
         return true;
      return false;
   }

   template<class T>
   struct Solver {
      using equation = map<int, T>;

      bool EquationSub(equation& a, const equation& e, T coeff) {
         for (const auto& q : e) {
            T v = a[q.first];
            T sub = coeff * q.second;
            if (ceq(v, sub))
               if (q.first)
                  a.erase(q.first);
               else
                  a[0] = 0;
            else a[q.first] -= coeff * q.second;
         }
         a[0] = a[0];
         if (a.size() < 3) {
            if (a.size() == 1) {
               if (abs(a[0]) > Eps<T>::eps) {
                  cerr << "Equation cannot be satisfied, aborted." << endl;
                  abort();
               }
               else return true;
            }
            if (a.size() == 2) {
               auto q = a.lower_bound(1);
               int id = q->first;
               T val = q->second;
               if (!solved.count(id)) {
                  T qval = -a[0] / val;
                  solved[id] = qval;
               }
               return true;
            }
         }
         return false;
      }

      inline static int varCounter = 0;
      using solver_t = Solver<T>;
      inline static solver_t* solver = nullptr;
      list<equation> equations;

      static int getVar() { return getSolver()->getNewVar(); }

      vector<int> relationMap;
      int getNewVar() {
         int n = ++varCounter;
         relationMap.push_back(n);
         return n;
      }
      int find(int n) { return relationMap[n] == n ? n : (relationMap[n] = find(relationMap[n])); }
      void merge(int u, int v) {
         relationMap[find(v)] = find(u);
      }

      static solver_t* getSolver() {
         return solver = solver ? solver : (new solver_t);
      }
      map<int, T> solved;
      T getSolution(int id) {
         if (solved.count(id))
            return solved[id];
         int u = find(id);
         using eqPtr = typename list<equation>::iterator;
         vector<eqPtr> involvedEquations;
         vector<int> eqErase;
         set<int> involvedIndices;

         for (auto l = equations.begin(), r = equations.end(); l != r; l++) {
            for (const auto& t : *l) {
               if (t.first && find(t.first) == u) {
                  involvedEquations.push_back(l);
                  eqErase.push_back(0);
                  break;
               }
            }
         }

         for (auto ep : involvedEquations) {
            for (const auto& t : *ep)
               if (t.first) involvedIndices.insert(t.first);
         }

         for (int ind : involvedIndices) {
            for (int i = 0, _ = involvedEquations.size(); i < _; ++i) {
               if (involvedEquations[i]->count(ind) && involvedEquations[i]->at(ind) != 0) {
                  if (i) {
                     swap(involvedEquations[0], involvedEquations[i]);
                     swap(eqErase[0], eqErase[i]);
                     break;
                  }
               }
            }
            if (involvedEquations[0]->count(ind) && involvedEquations[0]->at(ind) != 0) {
               T cc = involvedEquations[0]->at(ind);
               for (int i = 1, _ = involvedEquations.size(); i < _; ++i) {
                  T q = involvedEquations[i]->count(ind) ? involvedEquations[i]->at(ind) : 0;
                  if (q) {
                     bool erase = EquationSub(*(involvedEquations[i]), *(involvedEquations[0]), q / cc);
                     eqErase[i] = eqErase[i] || erase;
                  }
               }
            }
         }

         for (int i = 0, _ = involvedEquations.size(); i < _; ++i) {
            if (eqErase[i])
               equations.erase(involvedEquations[i]);
         }

         if (solved.count(id))
            return solved[id];
         else {
            cerr << "Variable cannot be solved, aborted. " << endl;
            abort();
         }
      }
      Solver(){
         relationMap.push_back(0);
      }
      void AddEquation(const equation& eq) {
         equation r;
         T coeffR = 0;
         int qid = 0;
         for (const auto& c : eq)
            if (c.first != 0) {
               if(solved.count(c.first)) coeffR += solved[c.first] * c.second;
               else {
                  r[c.first] = c.second;
                  if (qid) merge(qid, c.first);
                  else qid = c.first;
               }
            }else{
               coeffR -= c.second;
            }
         if (r.size() == 0 && abs(coeffR) > Eps<T>::eps) {
            cerr << "Equation cannot be satisfied, aborted." << endl;
            abort();
         }
         else {
            r[0] -= coeffR;
            equations.emplace_back(r);
         }
      }
      void printVar(int id) {
         if (id > 26) cout << "Z_" << id;
         else cout << names[id - 1];
      }
      bool printTerm(int id, T coeff, bool first) {
         if (coeff == 0) return first;
         if (coeff > 0) {
            if (!first) cout << " + ";
         }
         else {
            cout << " - ";
            coeff = -coeff;
         }
         if (coeff != 1) cout << coeff << " ";
         printVar(id);
         return false;
      }
      void PrintEquations() {
         for (const auto& v : solved) {
            T q = -v.second;
            if (q == 0) q = 0;
            printVar(v.first);
            cout << " = " << q << "\n";
         }
         for (const auto& eq : equations) {
            bool first = true;
            for (const auto& term : eq)
               if (term.first != 0)
                  first = printTerm(term.first, term.second, first);
            T q = 0;
            if (eq.count(0)) {
               q = -eq.at(0);
               if (q == 0) q = 0;
            }
            cout << " = " << q << "\n";
         }
      }
      static void addEq(const equation& eq) { solver_t::getSolver()->AddEquation(eq); }
      struct Expr {
         map<int, T> expr;
         Expr() {
            int id = solver_t::getVar();
            expr[id] = 1;
         }
         Expr(const map<int, T>& cp) {
            expr = cp;
         }
         Expr(T p) {
            expr[0] = p;
         }
         Expr(const Expr& c) {
            expr = c.expr;
         }
         Expr operator*(T b)const{
            Expr ret(expr);
            for (auto& p : ret)
               p.second *= b;
            return ret;
         }
         friend Expr operator*(T b, const Expr&a){
            Expr ret(a.expr);
            for (auto& p : ret)
               p.second *= b;
            return ret;
         }
         Expr operator+(T b)const {
            Expr ret(expr);
            ret.expr[0] += b;
            return ret;
         }
         friend Expr operator+(T b,const Expr&a){
            Expr ret(a.expr);
            ret.expr[0] += b;
            return ret;
         }

         Expr operator-(T b)const {
            Expr ret(expr);
            ret.expr[0] -= b;
            return ret;
         }
         friend Expr operator-(T b, const Expr& a){
            Expr ret(a.expr);
            for (auto& p : ret)
               p.second = -p.second;
            ret.expr[0] += b;
            return ret;
         }
         Expr operator-()const {
            Expr ret(expr);
            for (auto& p : ret)
               p.second = -p.second;
            return ret;
         }
         Expr operator+(const Expr& t)const {
            Expr ret(expr);
            for (const auto& i : t.expr)
               ret.expr[i.first] += i.second;
            return ret;
         }
         Expr operator-(const Expr& t)const {
            Expr ret(expr);
            for (const auto& i : t.expr)
               ret.expr[i.first] -= i.second;
            return ret;
         }
         Expr operator=(const Expr& t)const {
            Expr cret = (*this) - t;
            solver_t::addEq(cret.expr);
            return cret;
         }
         Expr operator=(T t)const {
            Expr cret = (*this) - t;
            solver_t::addEq(cret.expr);
            return cret;
         }
         operator T()const {
            T res = 0;
            for (const auto& i : expr) {
               if (i.first) res += solver_t::getSolver()->getSolution(i.first) * i.second;
               else res += i.second;
            }
            return res;
         }
         T operator*(const Expr& t) {
            return T(t) * T(*this);
         }
      };
   };
   using rsolve = Solver<double>;
   using var = rsolve::Expr;
}
