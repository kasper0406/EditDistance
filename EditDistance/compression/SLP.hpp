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

using namespace std;

namespace Compression {
  namespace SLP {
    enum ProductionType { KEY, LEFT, RIGHT };
    
    class Visitor {
    public:
      virtual void visit(class Terminal* terminal) = 0;
      virtual void visit(class NonTerminal* nonTerminal) = 0;
    };
    
    class Production {
    public:
      Production() : associatedString(""), derivedStringLength(0), isAddedInPartition(false), DISTTableIndex(-1), TYPE2_next(nullptr), reclaimCount_(0) { }
      
      virtual ~Production() { };
      virtual void accept(Visitor* visitor) = 0;
      
      /**
       * Information used for x-partition.
       * Consider refactoring this!
       */
      string associatedString;
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
      
    private:
      int64_t reclaimCount_;
    };
    
    class Terminal : public Production {
    public:
      Terminal(char symbol) : symbol_(symbol) { }
      
      void accept(Visitor* visitor) {
        visitor->visit(this);
      }
      
      char symbol() const { return symbol_; }
      
    private:
      char symbol_;
    };
    
    class NonTerminal : public Production {
    public:
      NonTerminal(int64_t name, Production* left, Production* right)
      : left_(left), right_(right), name_(name)
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
      }
      
      void setRight(Production* production) {
        if (right_ != nullptr)
          right_->decReclaimCount();
        
        production->incReclaimCount();
        right_ = production;
      }
      
      int64_t name() const { return name_; }
      
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
      
    private:
      Production* root_;
      int64_t derivedLength_; // Consider computing this value instead
      int64_t production_count_;
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
        
        ss_ << "\"n" << nonTerminal << "\" [label=\"X" << nonTerminal->name() << "\"];" << endl;
        ss_ << "\"n" << nonTerminal << "\" -- \"n" << nonTerminal->left() << "\";" << endl;
        nonTerminal->left()->accept(this);
        ss_ << "\"n" << nonTerminal << "\" -- \"n" << nonTerminal->right() << "\";" << endl;
        nonTerminal->right()->accept(this);
      }
      
      static string toDot(const SLP& slp) {
        SLPFoldedPrinter printer(slp);
        
        printer.ss_ << "graph G {" << endl;
        printer.slp_.root()->accept(&printer);
        printer.ss_ << "}" << endl;
        
        return printer.ss_.str();
      }
      
    private:
      SLPFoldedPrinter(const SLP& slp) : slp_(slp) { }
      
      const SLP& slp_;
      stringstream ss_;
      unordered_set<Production*> seen;
    };
    
    class Partitioner : public Visitor {
    public:
      void visit(Terminal* terminal) {
        assert(terminal->associatedString.empty() || terminal->associatedString == string(1, terminal->symbol()));
        terminal->derivedStringLength = 1;
        terminal->associatedString = string(1, terminal->symbol());
        
        if (x_ == 1) {
          addKeyProduction(terminal);
        }
      }
      
      void visit(NonTerminal* nonTerminal) {
        if (nonTerminal->derivedStringLength != 0) {
          if (nonTerminal->derivedStringLength >= x_ &&
              nonTerminal->left()->derivedStringLength < x_ && nonTerminal->right()->derivedStringLength < x_) {
            assert(nonTerminal->derivedStringLength < 2 * x_);
            
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
          stringstream ss;
          ss << nonTerminal->left()->associatedString << nonTerminal->right()->associatedString;
          
          assert(nonTerminal->associatedString.empty() || nonTerminal->associatedString == ss.str());
          nonTerminal->associatedString = ss.str();
          
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
          string S = "";
          Production* prev_selected = nullptr;
          
          while (S.size() < x_ && cur != stop) {
            if (cur->right() != prev) {
              assert(!cur->right()->associatedString.empty());
              assert(cur->right()->associatedString.length() < x_);
              
              stringstream ss;
              ss << S << cur->right()->associatedString;
              S = ss.str();
              
              cur->TYPE2_next = prev_selected;
              prev_selected = cur;
            }
            
            prev = cur;
            if (parentIterator != stack.rend()) {
              cur = parentIterator->second; ++parentIterator;
            } else {
              cur = nullptr;
            }
          }
          
          if (!S.empty()) {
            assert(prev_selected->associatedString.empty() || prev_selected->associatedString == S);
            prev_selected->associatedString = S;
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
          string S = "";
          Production* prev_selected = nullptr;
          
          while (S.size() < x_ && cur != stop) {
            if (cur->left() != prev) {
              assert(!cur->left()->associatedString.empty());
              assert(cur->left()->associatedString.length() < x_);
              
              stringstream ss;
              ss << cur->left()->associatedString << S;
              S = ss.str();
              
              cur->TYPE2_next = prev_selected;
              prev_selected = cur;
            }
            
            prev = cur;
            if (parentIterator != stack.rend()) {
              cur = parentIterator->second; ++parentIterator;
            } else {
              cur = nullptr;
            }
          }
          
          if (!S.empty()) {
            assert(prev_selected->associatedString.empty() || prev_selected->associatedString == S);
            prev_selected->associatedString = S;
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
        partition_.push_back(production);
        if (!production->isAddedInPartition) {
          if (type == LEFT)
            assert(!((NonTerminal*)production)->right()->associatedString.empty());
          
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
