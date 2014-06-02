#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include <list>
#include <tuple>

#include <sdsl/suffix_trees.hpp>

namespace Compression {
  namespace SLP {
    using namespace std;
    using namespace sdsl;
    
    enum ProductionType { KEY, LEFT, RIGHT };
    
    class Visitor {
    public:
      virtual void visit(class Terminal* terminal) = 0;
      virtual void visit(class NonTerminal* nonTerminal) = 0;
    };
    
    class Production {
    public:
      Production() : Production(0) { }
      
      Production(int64_t height) : Production(0, 1) { }
      
      Production(int64_t height, int64_t lzlen) : associatedStringLen(0), derivedStringLength(0), isAddedInPartition(false), DISTTableIndex(-1), TYPE2_next(nullptr), lzbuilder_length(lzlen), reclaimCount_(0), height_(height) { }
      
      virtual ~Production() { };
      virtual void accept(Visitor* visitor) = 0;
      
      /**
       * Information used for x-partition.
       * Consider refactoring this!
       */
      // string associatedString;
      uint64_t associatedStringLen;
      int64_t derivedStringLength;
      
      bool isAddedInPartition;
      int64_t DISTTableIndex; // Index into the dist table, making constant lookup possible
      
      ProductionType type;
      /////////////////////////////////////////
      // TYPE 2 information
      Production* TYPE2_next;
      /////////////////////////////////////////
      
      void incReclaimCount() { reclaimCount_++; }
      bool decReclaimCount() { return --reclaimCount_ <= 0; }
      int64_t reclaimCount() const { return reclaimCount_; }
      
      int64_t height() const {
        return height_;
      }
      
      void setHeight(int64_t height) {
        height_ = height;
      }
      
      virtual int8_t balance() const = 0;
      
      int64_t lzbuilder_length;
      
    private:
      int64_t reclaimCount_;
      
      int64_t height_;
    };
    
    class Terminal : public Production {
    public:
      Terminal(char symbol) : symbol_(symbol) { }
      
      void accept(Visitor* visitor) {
        visitor->visit(this);
      }
      
      char symbol() const { return symbol_; }
      
      int8_t balance() const {
        return 0;
      }
      
    private:
      char symbol_;
    };
    
    class NonTerminal : public Production {
    public:
      NonTerminal(int64_t name, Production* left, Production* right)
        : NonTerminal(name, left, right, max(left->height(), right->height()) + 1, left->lzbuilder_length + right->lzbuilder_length)
      { }
      
      NonTerminal(int64_t name, Production* left, Production* right, int64_t height, int64_t lzlen)
        : Production(height, lzlen), left_(left), right_(right), name_(name)
      {
        left_->incReclaimCount();
        right_->incReclaimCount();
      }
      
      ~NonTerminal() {
        if (left_->decReclaimCount())
          delete left_;
        if (right_->decReclaimCount())
          delete right_;
        left_ = right_ = nullptr;
      }
      
      NonTerminal(NonTerminal&& other) {
        left_ = move(other.left_);
        right_ = move(other.right_);
        name_ = move(other.name_);
      }
      
      NonTerminal(const NonTerminal& other) {
        throw runtime_error("Do not copy non-terminals!");
      }
      
      NonTerminal& operator=(const NonTerminal&& other) {
        left_ = move(other.left_);
        right_ = move(other.right_);
        name_= move(other.name_);
        return *this;
      }
      
      NonTerminal& operator=(const NonTerminal& other) {
        throw runtime_error("Do not copy non-terminals!");
      }
      
      void accept(Visitor* visitor) {
        visitor->visit(this);
      }
      
      Production* left() const { return left_; }
      Production* right() const { return right_; }
      
      void setLeft(Production* production) {
        if (left_ != nullptr)
          left_->decReclaimCount();
        
        production->incReclaimCount();
        left_ = production;
        
        update_height();
      }
      
      void setRight(Production* production) {
        if (right_ != nullptr)
          right_->decReclaimCount();
        
        production->incReclaimCount();
        right_ = production;
        
        update_height();
      }
      
      int64_t name() const { return name_; }
      void setName(int64_t name) { name_ = name; }
      
      int8_t balance() const {
        return left_->height() - right_->height();
      }
      
      void update_height() {
        this->setHeight(max(left_->height(), right_->height()) + 1);
        this->lzbuilder_length = left_->lzbuilder_length + right_->lzbuilder_length;
      }
      
    private:
      Production *left_, *right_;
      int64_t name_;
    };
    
    class SLP {
    public:
      SLP(Production* root, int64_t derivedLength, int64_t production_count)
        : root_(root), derivedLength_(derivedLength), production_count_(production_count)
      { }
      
      ~SLP() {
        if (root_ != nullptr)
          delete root_;
      }
      
      SLP(SLP&& other) {
        root_ = move(other.root_);
      }
      
      SLP& operator=(SLP&& other) {
        root_ = move(other.root_);
        return *this;
      }
      
      SLP(const SLP& other) {
        throw runtime_error("Do not copy SLPs!");
      }
      
      SLP& operator=(const SLP& other) {
        throw runtime_error("Do not copy SLPs!");
      }
      
      Production* root() const {
        return root_;
      }
      
      void setRoot(Production* root) {
        root->decReclaimCount();
        root_->incReclaimCount();
        
        root_ = root;
      }
      
      int64_t productions() const {
        return production_count_;
      }
      
      void setProductions(int64_t count) {
        production_count_ = count;
      }
      
      int64_t derivedLength() const {
        return derivedLength_;
      }
      
      void setDerivedLength(int64_t length) {
        derivedLength_ = length;
      }
      
      int64_t height() const {
        return root_->height();
      }
      
      double compressionFactor() const {
        return double(productions()) / derivedLength();
      }
      
    private:
      Production* root_;
      int64_t derivedLength_; // Consider computing this value instead
      int64_t production_count_;
    };
    
    class SLPUnfoldedPrinter : public Visitor {
    public:
      void visit(Terminal* terminal) {
        ss_ << "\"n" << ++count_ << "\" [label=\"" << terminal->symbol() << "\",shape=\"plaintext\"];" << endl;
      }
      
      void visit(NonTerminal* nonTerminal) {
        int64_t node_count = ++count_;
        ss_ << "\"n" << node_count << "\" [label=\"X" << nonTerminal->name() << "\"];" << endl;
        ss_ << "\"n" << node_count << "\" -- \"n" << node_count + 1 << "\";" << endl;
        nonTerminal->left()->accept(this);
        ss_ << "\"n" << node_count << "\" -- \"n" << count_ + 1 << "\";" << endl;
        nonTerminal->right()->accept(this);
      }
      
      static string toDot(const SLP& slp) {
        SLPUnfoldedPrinter printer(slp);
        
        printer.ss_ << "graph G {" << endl;
        printer.ss_ << "graph [ordering=\"out\"];" << endl;
        printer.slp_.root()->accept(&printer);
        printer.ss_ << "}" << endl;
        
        return printer.ss_.str();
      }
      
    private:
      SLPUnfoldedPrinter(const SLP& slp) : slp_(slp), count_(0) { }
      
      const SLP& slp_;
      stringstream ss_;
      int64_t count_;
    };
    
    class SLPFoldedPrinter : public Visitor {
    public:
      void visit(Terminal* terminal) {
        if (seen.count(terminal) > 0) return;
        seen.insert(terminal);
        
        ss_ << "\"n" << terminal << "\" [label=\"" << terminal->symbol() << "\",shape=\"plaintext\"];" << endl;
      }
      
      void visit(NonTerminal* nonTerminal) {
        if (seen.count(nonTerminal) > 0) return;
        seen.insert(nonTerminal);
        
        // ss_ << "\"n" << nonTerminal << "\" [label=\"X" << nonTerminal->name() << ", height=" << nonTerminal->height() << ", rc=" << nonTerminal->reclaimCount() << ", lzlen=" << nonTerminal->lzbuilder_length << "\"];" << endl;
        
        ss_ << "\"n" << nonTerminal << "\" [label=\"X" << nonTerminal->name() << "\"];" << endl;
        
        ss_ << "\"n" << nonTerminal << "\" -> \"n" << nonTerminal->left() << "\";" << endl;
        nonTerminal->left()->accept(this);
        ss_ << "\"n" << nonTerminal << "\" -> \"n" << nonTerminal->right() << "\";" << endl;
        nonTerminal->right()->accept(this);
      }
      
      static string toDot(const SLP& slp) {
        SLPFoldedPrinter printer;
        
        printer.ss_ << "digraph G {" << endl;
        printer.ss_ << "graph [ordering=\"out\"];" << endl;
        slp.root()->accept(&printer);
        printer.ss_ << "}" << endl;
        
        return printer.ss_.str();
      }
      
      static string toDot(Production* root) {
        SLPFoldedPrinter printer;
        
        printer.ss_ << "graph G {" << endl;
        root->accept(&printer);
        printer.ss_ << "}" << endl;
        
        return printer.ss_.str();
      }
      
    private:
      SLPFoldedPrinter() { }
      
      stringstream ss_;
      unordered_set<Production*> seen;
    };
    
    class SimpleSLPBuilder {
    public:
      static unique_ptr<SLP> build(string input) {
        SimpleSLPBuilder builder;
        auto root = builder.buildTree(input);
        return unique_ptr<SLP>(new SLP(root, input.size(), builder.production_count));
      }
      
      static string name() {
        return "TrivialSLP";
      }
      
    private:
      SimpleSLPBuilder() : production_count(0), count(0) { }
      
      Production* buildTree(string input) {
        production_count++;
        if (input.size() == 1) {
          return new Terminal(input[0]);
        } else {
          int64_t len = input.length();
          
          return new NonTerminal(++count,
                                 buildTree(input.substr(0, (len + 1) / 2)),
                                 buildTree(input.substr((len + 1) / 2, len / 2)));
        }
      }
      
      int64_t production_count;
      int64_t count;
    };
    
    class SimpleCompressionSLPBuilder {
    public:
      static unique_ptr<SLP> build(string input) {
        SimpleCompressionSLPBuilder builder;
        auto root = builder.buildTree(input);
        return unique_ptr<SLP>(new SLP(root, input.size(), builder.production_count));
      }
      
      static string name() {
        return "SimpleCompressionSLP";
      }
      
    private:
      SimpleCompressionSLPBuilder() : production_count(0), count(0) { }
      
      Production* buildTree(string input) {
        // Check if string is already derived by some variable
        if (compressed.count(input) > 0) {
          // Already in map
          return compressed[input];
        } else {
          Production* result;
          production_count++;
          
          if (input.size() == 1) {
            result = new Terminal(input[0]);
          } else {
            int64_t len = input.length();
            
            NonTerminal* nonTerminal = new NonTerminal(++count,
                                                       buildTree(input.substr(0, (len + 1) / 2)),
                                                       buildTree(input.substr((len + 1) / 2, len / 2)));
            
            result = nonTerminal;
          }
          
          compressed[input] = result;
          return result;
        }
      }
      
      int64_t production_count;
      int64_t count;
      unordered_map<string, Production*> compressed;
    };
    
    class StringDeriver : private Visitor {
    public:
      static string getDerivedString(Production* production) {
        StringDeriver deriver;
        production->accept(&deriver);
        return deriver.ss.str();
      }
      
    private:
      void visit(Terminal* t) {
        ss << t->symbol();
      }
      
      void visit(NonTerminal* nt) {
        nt->left()->accept(this);
        nt->right()->accept(this);
      }
      
      stringstream ss;
    };
    
    class LZFactorize {
    public:
      static vector<pair<uint64_t,uint64_t>> naive_lz_factorize(string str) {
        // cout << str << endl;
        
        vector<pair<uint64_t,uint64_t>> factors;
        
        uint64_t i = 0;
        while (i < str.length()) {
          int64_t length = 1;
          uint64_t found_at = i;
          for (; length <= str.length() - i; ++length) {
            uint64_t found = str.substr(0, i).find(str.substr(i, length));
            if (found != string::npos) {
              found_at = found;
            } else {
              break;
            }
          }
          
          uint64_t begin = found_at;
          uint64_t end = found_at + (max(length - 2, int64_t(0)));
          
          factors.push_back({ begin, end });
          // cout << str.substr(begin, end - begin + 1) << endl;
          
          i += end - begin + 1;
        }
        
        return factors;
      }
            
      static vector<pair<uint64_t,uint64_t>> lz_factorize(string str) {
        cst_sct3<> cst;
        construct_im(cst, str.c_str(), 1);
        vector<pair<uint64_t,uint64_t>> factors;
        
        vector<uint64_t> first_time_seen(cst.nodes(), numeric_limits<uint64_t>::max());
        for (auto iter = cst.begin(); iter != cst.end(); ++iter) {
          auto current_id = cst.id(*iter);
          
          // cout << current_id << endl;
          
          if (cst.is_leaf(*iter)) {
            first_time_seen[current_id] = cst.sn(*iter); // str.length() - cst.depth(*iter);
          }
          
          if (cst.root() != *iter && first_time_seen[current_id] != numeric_limits<uint64_t>::max()) {
            auto parent_id = cst.id(cst.parent(*iter));
            first_time_seen[parent_id] = min(first_time_seen[parent_id], first_time_seen[current_id]);
          }
        }
        
        for (uint64_t factor_start = 0; factor_start < str.length();) {
          auto current_node = cst.root(), previous_node = cst.root();
          uint64_t current_edge_length = 0, current_edge_matched = 0;
          
          for (uint64_t factor_len = 1; ; ++factor_len) {
            if (current_edge_matched >= current_edge_length) {
              previous_node = current_node;
              current_node = cst.child(current_node, str[factor_start + factor_len - 1]);
              current_edge_length = cst.depth(current_node) - cst.depth(previous_node);
              current_edge_matched = 0;
            }
            
            if (first_time_seen[cst.id(current_node)] + factor_len - 1 >= factor_start) {
              if (factor_len == 1) {
                factors.push_back({ factor_start, factor_start });
                factor_start++;
              } else {
                const auto pos = first_time_seen[cst.id(previous_node)];
                factors.push_back({ pos, factor_len - 2 + pos });
                factor_start += factor_len - 1;
              }
              
              break;
            }
            
            previous_node = current_node;
            current_edge_matched++;
          }
        }
        
        return factors;
      }
    };
    
    class LZSLPBuilder {
    public:
      static unique_ptr<SLP> build(string input) {
        Production* root = nullptr;
        
        // cout << "Factors:" << endl;
        
        int64_t processed = 0;
        auto factors = LZFactorize::lz_factorize(input);
        for (auto factor : factors) {
          // cout << input.substr(get<0>(factor), get<1>(factor) - get<0>(factor) + 1) << endl;
          
          if (get<0>(factor) == processed) {
            assert(get<1>(factor) == get<0>(factor));
            Terminal* term = new Terminal(input[get<0>(factor)]);
            if (root != nullptr) {
              root = concat(root, term);
            } else {
              root = term;
            }
          } else {
            auto productions = NodeHarvester::productions(root, get<0>(factor), get<1>(factor));
            
            vector<Production*> left, right;
            left.reserve(productions.size()); right.reserve(productions.size());
            int64_t prev_height = 0;
            for (auto production : productions) {
              if (production->height() >= prev_height) {
                assert(right.empty());
                left.push_back(production);
                prev_height = production->height();
              } else {
                right.push_back(production);
              }
            }
            
            // cout << "ls = " << left.size() << endl;
            // cout << "rs = " << right.size() << endl;
            
            root = concat(root, concat(concat(left), concat(right)));
          }
          
          processed += get<1>(factor) - get<0>(factor) + 1;
        }
        
        int64_t num_productions = CounterAndNamer::run(root);
        assert(input == StringDeriver::getDerivedString(root));
        
        return unique_ptr<SLP>(new SLP(root, input.size(), num_productions));
      }
      
      static string name() {
        return "LZSLPBuilder";
      }
      
    private:
      static Production* concat(Production* A, Production* B) {
        if (A == nullptr) return B;
        else if (B == nullptr) return A;
        
//        cout << "A height: " << A->height() << endl;
//        cout << "B height: " << B->height() << endl;
        
        Production *root, *low_tree;
        bool A_is_root = A->height() >= B->height();
        if (A_is_root) {
          root = A;
          low_tree = B;
        } else {
          root = B;
          low_tree = A;
        }
        
        vector<NonTerminal*> parents;
        Production *v = root;
        while (v->height() > low_tree->height() + 1) {
          auto v_nonterm = static_cast<NonTerminal*>(v);
          
          assert(v->reclaimCount() >= 1 || root == v);
          if (root != v) {
            assert(!parents.empty());
            
            // Other guys are pointing to this ancestor. We must not alter the subtree!
            // Solution: Create a new node.
            NonTerminal* v_ = new NonTerminal(0, v_nonterm->left(), v_nonterm->right());
            if (A_is_root) {
              parents.back()->setRight(v_);
            } else {
              parents.back()->setLeft(v_);
            }
            
            v = v_nonterm = v_;
          } else {
            root = v = v_nonterm = new NonTerminal(0, v_nonterm->left(), v_nonterm->right());
          }
          
          // Static cast is safe, since terminals will alway have height 0.
          parents.push_back(v_nonterm);
          if (A_is_root) {
            v = v_nonterm->right();
          } else {
            v = v_nonterm->left();
          }
        }
        
        if (!parents.empty()) {
          NonTerminal* v_;
          if (A_is_root) {
            v_ = new NonTerminal(0, v, low_tree);
            parents.back()->setRight(v_);
          } else {
            v_ = new NonTerminal(0, low_tree, v);
            parents.back()->setLeft(v_);
          }
          assert(v_->balance() >= -1 && v_->balance() <= 1);
          
          // Update the height of parent path and do rotations if required
          for (auto iter = parents.rbegin(); iter != parents.rend(); ++iter) {
            auto ancestor = *iter;
            ancestor->update_height();
            
            // Check if ancestor is balanced.
//            cout << "Balance: " << int64_t(ancestor->balance()) << endl;
            assert(ancestor->balance() >= -2 && ancestor->balance() <= 2);
            if (ancestor->balance() == -2) {
              assert(A_is_root);
//              cout << "Rebalance required: Right too high." << endl;
              
              auto A = *iter;
              auto B = A->left();
              auto D = static_cast<NonTerminal*>(A->right()); // Cast alright, since height >= 2
              auto C = D->left();
              auto E = D->right();
              
              NonTerminal* A_;
              if (E->height() > C->height()) {
                // Rotation 1
                auto D_ = new NonTerminal(0, B, C);
                A_ = new NonTerminal(0, D_, E);
              } else {
                // Rotation 2
                auto F = static_cast<NonTerminal*>(C)->left();
                auto G = static_cast<NonTerminal*>(C)->right();
                
                auto C_ = new NonTerminal(0, B, F);
                auto D_ = new NonTerminal(0, G, E);
                A_ = new NonTerminal(0, C_, D_);
              }
              
              auto next_ancestor_iter = next(iter);
              if (next_ancestor_iter == parents.rend()) {
                return A_;
              } else {
                auto next_ancestor = *next_ancestor_iter;
                assert(next_ancestor->right() == A);
                next_ancestor->setRight(A_);
              }
            } else if (ancestor->balance() == 2) {
              assert(!A_is_root);
//              cout << "Rebalance required: Left too high." << endl;
              
              auto A = *iter;
              auto B = static_cast<NonTerminal*>(A->left()); // Cast alright, due to balance
              auto E = A->right();
              auto C = B->left();
              auto D = B->right();
              
              NonTerminal* A_;
              if (C->height() > D->height()) {
                // Rotation 1
                auto B_ = new NonTerminal(0, D, E);
                A_ = new NonTerminal(0, C, B_);
              } else {
                // Rotation 2
                auto F = static_cast<NonTerminal*>(D)->left();
                auto G = static_cast<NonTerminal*>(D)->right();
                
                auto B_ = new NonTerminal(0, C, F);
                auto D_ = new NonTerminal(0, G, E);
                A_ = new NonTerminal(0, B_, D_);
              }
              
              auto next_ancestor_iter = next(iter);
              if (next_ancestor_iter == parents.rend()) {
                return A_;
              } else {
                auto next_ancestor = *next_ancestor_iter;
                assert(next_ancestor->left() == A);
                next_ancestor->setLeft(A_);
              }
            }
          }
          
          return static_cast<NonTerminal*>(root);
        } else {
          root = new NonTerminal(0, A, B);
//          cout << "Root balance: " << int64_t(root->balance()) << endl;
          assert(root->balance() >= -1 && root->balance() <= 1);
          return static_cast<NonTerminal*>(root);
        }
        
        throw runtime_error("Should not happen.");
      }
      
      static Production* concat(vector<Production*> productions) {
        if (productions.size() == 0) return nullptr;
        if (productions.size() == 1) return productions.front();
        
        Production* concatted = productions.front();
        for (auto iter = next(productions.begin()); iter != productions.end(); ++iter) {
          concatted = concat(concatted, *iter);
        }
        
        return concatted;
      }
      
      class NodeHarvester : private Visitor {
      public:
        static vector<Production*> productions(Production* root, int64_t start, int64_t end) {
          NodeHarvester harvester(start, end);
          root->accept(&harvester);
          
          return harvester.productions_;
        }
        
      private:
        NodeHarvester(int64_t start, int64_t end) : start_(start), end_(end) { }
        
        void visit(Terminal* terminal) {
          assert(start_ == end_);
          productions_.push_back(terminal);
        }
        
        void visit(NonTerminal* nonTerminal) {
          assert(start_ <= end_);
          
          if (nonTerminal->lzbuilder_length == end_ - start_ + 1) {
            productions_.push_back(nonTerminal);
            return;
          }
          
          const int64_t left_len = nonTerminal->left()->lzbuilder_length;
          if (left_len > start_) {
            // Call recursively on left
            int64_t tmp_end = end_;
            end_ = min(end_, left_len - 1);
            nonTerminal->left()->accept(this);
            end_ = tmp_end;
          }
          if (left_len <= end_) {
            // Call recursively on right
            end_ -= left_len;
            start_ = max((int64_t)0, start_ - left_len);
            nonTerminal->right()->accept(this);
            end_ += left_len;
            start_ += left_len;
          }
        }
        
        vector<Production*> productions_;
        int64_t start_, end_;
      };
      
      class CounterAndNamer : private Visitor {
      public:
        // Returns: #productions
        static int64_t run(Production* root) {
          CounterAndNamer can;
          root->accept(&can);
          return can.counter;
        }
        
      private:
        CounterAndNamer() : counter(0) { }
        
        void visit(Terminal* terminal) {
          if (seen.count(terminal) > 0) return;
          seen.insert(terminal);
          
          ++counter;
        }
        
        void visit(NonTerminal* nonTerminal) {
          if (seen.count(nonTerminal) > 0) return;
          seen.insert(nonTerminal);
          
          ++counter;
          nonTerminal->setName(counter);
          
          nonTerminal->left()->accept(this);
          nonTerminal->right()->accept(this);
        }
        
        unordered_set<Production*> seen;
        int64_t counter;
      };
    };
    
    /**
     * Blow up the string generated by the SLP by replacing all terminals with non-terminals producting original terminal and a '$' symbol.
     * This is used by the LCS EditDistance calculator.
     */
    class BlowUpSLPTransformer : private Visitor {
    public:
      static void blowUpSLP(SLP* slp) {
        BlowUpSLPTransformer transformer(slp);
        slp->root()->accept(&transformer);
        
        slp->setProductions(slp->productions() + transformer.transformed.size() + 1);
        slp->setDerivedLength(slp->derivedLength() * 2);
      }
      
    private:
      BlowUpSLPTransformer(SLP* slp) {
        slp_ = slp;
        count_ = slp->productions();
        special_symbol_ = new Terminal('$');
        cur_parent_ = nullptr;
      }
      
      void visit(Terminal* terminal) {
        if (seen.count(terminal) == 0) {
          seen.insert(terminal);
          
          NonTerminal* nonTerminal = new NonTerminal(++count_, terminal, special_symbol_);
          seen.insert(nonTerminal);
          transformed.insert({ terminal, nonTerminal });
        }
        
        if (cur_parent_ == nullptr) {
          slp_->setRoot(transformed[terminal]);
        } else if (cur_parent_->left() == terminal) {
          cur_parent_->setLeft(transformed[terminal]);
        } else {
          assert(cur_parent_->right() == terminal);
          cur_parent_->setRight(transformed[terminal]);
        }
      }
      
      void visit(NonTerminal* nonTerminal) {
        if (seen.count(nonTerminal) > 0) return;
        seen.insert(nonTerminal);
        
        cur_parent_ = nonTerminal;
        nonTerminal->left()->accept(this);
        
        cur_parent_ = nonTerminal;
        nonTerminal->right()->accept(this);
      }
      
      unordered_set<Production*> seen;
      unordered_map<Terminal*, NonTerminal*> transformed;
      
      int64_t count_;
      NonTerminal* cur_parent_;
      Terminal* special_symbol_;
      SLP* slp_;
    };
    
    class Partitioner : public Visitor {
    public:
      void visit(Terminal* terminal) {
        // assert(terminal->associatedString.empty() || terminal->associatedString == string(1, terminal->symbol()));
        terminal->derivedStringLength = 1;
        // terminal->associatedString = string(1, terminal->symbol());
        terminal->associatedStringLen = 1;
        
        if (x_ == 1) {
          addKeyProduction(terminal);
        }
      }
      
      void visit(NonTerminal* nonTerminal) {
        if (nonTerminal->derivedStringLength != 0) {
          if (nonTerminal->derivedStringLength >= x_ &&
              nonTerminal->left()->derivedStringLength < x_ && nonTerminal->right()->derivedStringLength < x_) {
            assert(nonTerminal->derivedStringLength < 2 * x_);
            assert(nonTerminal->associatedStringLen == nonTerminal->derivedStringLength);
            
            addKeyProduction(nonTerminal);
            return; // We are done with the subtree
          }
        }
        
        parentStack_.push_back({ LEFT, nonTerminal });
        nonTerminal->left()->accept(this);
        parentStack_.pop_back();
        
        parentStack_.push_back({ RIGHT, nonTerminal });
        nonTerminal->right()->accept(this);
        parentStack_.pop_back();
        
        const int64_t leftLen = nonTerminal->left()->derivedStringLength,
        rightLen = nonTerminal->right()->derivedStringLength;
        nonTerminal->derivedStringLength = leftLen + rightLen;
        if (leftLen < x_ && rightLen < x_) {
          // stringstream ss;
          // ss << nonTerminal->left()->associatedString << nonTerminal->right()->associatedString;
          
          // assert(nonTerminal->associatedString.empty() || nonTerminal->associatedString == ss.str());
          // nonTerminal->associatedString = ss.str();
          nonTerminal->associatedStringLen = nonTerminal->left()->associatedStringLen + nonTerminal->right()->associatedStringLen;
          
          if (leftLen + rightLen >= x_) {
            addKeyProduction(nonTerminal);
          }
        }
      }
      
      /**
       * Harvest strings in between key vertices i and i + 1
       */
      void harvestIntermediateNodes(Production* prevKeyVertex, Production* newKeyVertex) {
        // Find the LCA of i and i + 1
        int64_t k = 0;
        while (k < previousParentStack_.size() && k < parentStack_.size() && previousParentStack_[k] == parentStack_[k]) ++k;
        if (!(k < previousParentStack_.size() && k < parentStack_.size() && previousParentStack_[k].second == parentStack_[k].second))
          k--;
        Production* lca;
        if (k < 0) lca = slp_.root();
        else if (previousParentStack_.size() > k) lca = previousParentStack_[k].second;
        else if (parentStack_.size() > k) lca = parentStack_[k].second;
        else lca = slp_.root();
        
        harvestLeftPath(prevKeyVertex, previousParentStack_, lca);
        harvestRightPath(newKeyVertex, parentStack_, lca);
      }
      
      // Returns (partition, blocks)
      static pair<vector<Production*>, vector<Production*>> partition(const SLP& slp, int64_t x) {
        Partitioner partitioner(x, slp);
        
        // Find key vertices
        partitioner.partition_.reserve(slp.derivedLength() / x);
        partitioner.blocks_.reserve(slp.derivedLength() / x);
        
        slp.root()->accept(&partitioner);
        partitioner.harvestLeftPath(partitioner.partition_[partitioner.partition_.size() - 1],
                                    partitioner.previousParentStack_, nullptr); // Right of last key vertex
        
        return { partitioner.partition_, partitioner.blocks_ };
      }
      
    private:
      Partitioner(int64_t x, const SLP& slp) : x_(x), slp_(slp) { }
      
      void harvestLeftPath(Production* start, const vector<pair<ProductionType, NonTerminal*>>& stack, const Production* stop) {
        if (start == stop || stack.size() == 0) return;
        
        auto parentIterator = stack.rbegin();
        NonTerminal* cur = parentIterator->second; ++parentIterator;
        Production* prev = start;
        while (cur != stop) {
          // string S = "";
          int64_t S_len = 0;
          Production* prev_selected = nullptr;
          // Production* prev_added_to_partition = nullptr;
          
          while (S_len < x_ && cur != stop) {
            if (cur->right() != prev) {
              assert(cur->right()->associatedStringLen != 0);
              assert(cur->right()->associatedStringLen < x_);
              
              // S += cur->right()->associatedString;
              S_len += cur->right()->associatedStringLen;
              
              /*
              stringstream ss;
              ss << S << cur->right()->associatedString;
              S = ss.str();
               */
              
              cur->TYPE2_next = prev_selected;
              cur->type = LEFT;
              
              prev_selected = cur;
            }
            
            prev = cur;
            if (parentIterator != stack.rend()) {
              cur = parentIterator->second; ++parentIterator;
            } else {
              cur = nullptr;
            }
          }
          
          if (S_len != 0) {
            // assert(prev_selected->associatedString.empty() || prev_selected->associatedString == S);
            
            // prev_selected->TYPE2_next = prev_added_to_partition;
            // prev_added_to_partition = prev_selected;
            
            prev_selected->associatedStringLen = S_len;
            addToPartition(prev_selected, LEFT);
          }
        }
      }
      
      void harvestRightPath(Production* start, const vector<pair<ProductionType, NonTerminal*>>& stack, const Production* stop) {
        if (start == stop || stack.size() == 0) return;
        
        auto parentIterator = stack.rbegin();
        NonTerminal* cur = parentIterator->second; ++parentIterator;
        Production* prev = start;
        while (cur != stop) {
          // string S = "";
          int64_t S_len = 0;
          Production* prev_selected = nullptr;
          // Production* prev_added_to_partition = nullptr;
          
          while (S_len < x_ && cur != stop) {
            if (cur->left() != prev) {
              assert(cur->left()->associatedStringLen != 0);
              assert(cur->left()->associatedStringLen < x_);
              
              /*
              stringstream ss;
              ss << cur->left()->associatedString << S;
              S = ss.str();
               */
              S_len += cur->left()->associatedStringLen;
              
              cur->TYPE2_next = prev_selected;
              cur->type = RIGHT;
              
              prev_selected = cur;
            }
            
            prev = cur;
            if (parentIterator != stack.rend()) {
              cur = parentIterator->second; ++parentIterator;
            } else {
              cur = nullptr;
            }
          }
          
          if (S_len != 0) {
            // assert(prev_selected->associatedString.empty() || prev_selected->associatedString == S);
            
            // prev_selected->TYPE2_next = prev_added_to_partition;
            // prev_added_to_partition = prev_selected;
            
            // prev_selected->associatedString = S;
            prev_selected->associatedStringLen = S_len;
            addToPartition(prev_selected, RIGHT);
          }
        }
      }
      
      void addKeyProduction(Production* production) {
        if (partition_.size() == 0)
          harvestRightPath(production, parentStack_, nullptr);
        else {
          harvestIntermediateNodes(previousKeyVertex_, production);
        }
        
        addToPartition(production, KEY);
        
        previousParentStack_ = vector<pair<ProductionType, NonTerminal*>>(parentStack_);
        previousKeyVertex_ = production;
      }
      
      void addToPartition(Production* production, ProductionType type) {
        assert(production->associatedStringLen != 0);
        assert(production->associatedStringLen == production->derivedStringLength || type != KEY);
        
        partition_.push_back(production);
        if (!production->isAddedInPartition) {
          if (type == LEFT)
            assert(((NonTerminal*)production)->right()->associatedStringLen != 0);
          if (type == RIGHT)
            assert(((NonTerminal*)production)->left()->associatedStringLen != 0);
          
          production->isAddedInPartition = true;
          production->type = type;
          blocks_.push_back(production);
        }
      }
      
      const int64_t x_;
      const SLP& slp_;
      
      Production* previousKeyVertex_;
      vector<pair<ProductionType, NonTerminal*>> previousParentStack_;
      vector<pair<ProductionType, NonTerminal*>> parentStack_;
      
      vector<Production*> partition_;
      vector<Production*> blocks_;
    };
  }
}
