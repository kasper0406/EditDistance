#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

using namespace std;

namespace SLP {
  class Visitor {
  public:
    virtual void visit(class Terminal* terminal) = 0;
    virtual void visit(class NonTerminal* nonTerminal) = 0;
  };
  
  class Production {
  public:
    Production() : associatedString(""), derivedStringLength(0), parent_(nullptr) { }
    
    virtual ~Production() { };
    virtual void accept(Visitor* visitor) = 0;
    
    void setParent(NonTerminal* parent) { parent_ = parent; }
    NonTerminal* parent() const { return parent_; }
    
    
    /**
     * Information used for x-partition.
     * Consider refactoring this!
     */
    string associatedString;
    uint64_t derivedStringLength;
    
  protected:
    NonTerminal* parent_;
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
    NonTerminal(uint64_t name, Production* left, Production* right)
      : left_(left), right_(right), name_(name)
    { }
    
    ~NonTerminal() {
      delete left_;
      delete right_;
    }
    
    NonTerminal(NonTerminal&& other) {
      left_ = move(other.left_);
      right_ = move(other.right_);
      parent_ = move(other.parent_);
      name_= move(other.name_);
    }
    
    NonTerminal(const NonTerminal& other) {
      throw runtime_error("Do not copy non-terminals!");
    }
    
    NonTerminal& operator=(const NonTerminal&& other) {
      left_ = move(other.left_);
      right_ = move(other.right_);
      parent_ = move(other.parent_);
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
    
    uint64_t name() const { return name_; }
    
  private:
    Production *left_, *right_;
    uint64_t name_;
  };
  
  class SLP {
  public:
    SLP(Production* root, uint64_t derivedLength) : root_(root), derivedLength_(derivedLength) { }
    
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
    
    uint64_t derivedLength() const {
      return derivedLength_;
    }
    
  private:
    Production* root_;
    uint64_t derivedLength_; // Consider computing this value instead
  };
  
  class SimpleSLPBuilder {
  public:
    static SLP build(string input) {
      return SLP(SimpleSLPBuilder().buildTree(input), input.size());
    }
    
  private:
    SimpleSLPBuilder() : count(0) { }
    
    Production* buildTree(string input) {
      if (input.size() == 1) {
        return new Terminal(input[0]);
      } else {
        uint64_t len = input.length();
        
        NonTerminal* nonTerminal = new NonTerminal(++count,
                                                   buildTree(input.substr(0, (len + 1) / 2)),
                                                   buildTree(input.substr((len + 1) / 2, len / 2)));
        nonTerminal->right()->setParent(nonTerminal);
        nonTerminal->left()->setParent(nonTerminal);
        
        return nonTerminal;
      }
    }
    
    uint64_t count;
  };
  
  class SLPUnfoldedPrinter : public Visitor {
  public:
    void visit(Terminal* terminal) {
      ss_ << "\"n" << ++count_ << "\" [label=\"" << terminal->symbol() << "\",shape=\"plaintext\"];" << endl;
    }
    
    void visit(NonTerminal* nonTerminal) {
      uint64_t node_count = ++count_;
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
    uint64_t count_;
  };
  
  class Partitioner : public Visitor {
  public:
    void visit(Terminal* terminal) {
      terminal->derivedStringLength = 1;
      terminal->associatedString = string(1, terminal->symbol());
      
      if (x_ == 1) {
        keyVertices.push_back(terminal);
      }
    }
    
    void visit(NonTerminal* nonTerminal) {
      if (nonTerminal->derivedStringLength != 0) {
        assert(nonTerminal->derivedStringLength >= x_ && nonTerminal->derivedStringLength < 2 * x_);
        assert(nonTerminal->left()->derivedStringLength < x_ && nonTerminal->right()->derivedStringLength < x_);
        
        keyVertices.push_back(nonTerminal);
        
        return; // This production is already handled.
      }
      
      nonTerminal->left()->accept(this);
      nonTerminal->right()->accept(this);
      
      const uint64_t leftLen = nonTerminal->left()->derivedStringLength,
                     rightLen = nonTerminal->right()->derivedStringLength;
      nonTerminal->derivedStringLength = leftLen + rightLen;
      if (leftLen < x_ && rightLen < x_) {
        stringstream ss;
        ss << nonTerminal->left()->associatedString << nonTerminal->right()->associatedString;
        nonTerminal->associatedString = ss.str();
        
        if (leftLen + rightLen >= x_) {
          keyVertices.push_back(nonTerminal);
        }
      }
    }
    
    static NonTerminal* lca(Production* a, Production* b) {
      // TODO: Consider faster implementation
      assert(a != b);
      
      unordered_set<Production*> seen;
      seen.insert(a); seen.insert(b);
      
      NonTerminal* path1 = a->parent();
      NonTerminal* path2 = b->parent();
      while (path1 != nullptr || path2 != nullptr) {
        if (path1 != nullptr) {
          if (seen.count(path1) > 0) return path1;
          seen.insert(path1); path1 = path1->parent();
        }
        
        if (path2 != nullptr) {
          if (seen.count(path2) > 0) return path2;
          seen.insert(path2); path2 = path2->parent();
        }
      }
      
      throw runtime_error("This should not happen! Could not find LCA!");
    }
    
    void harvestLeftPath(Production* start, const Production* stop, vector<pair<string, Production*>>& partition) {
      NonTerminal* cur = start->parent();
      Production* prev = start;
      while (cur != stop) {
        string S = "";
        while (S.size() < x_ && cur != stop) {
          if (cur->right() != prev) {
            assert(!cur->right()->associatedString.empty());
            assert(cur->right()->associatedString.length() < x_);
            
            stringstream ss;
            ss << S << cur->right()->associatedString;
            S = ss.str();
          }
          
          prev = cur;
          cur = cur->parent();
        }
        
        if (!S.empty()) {
          partition.push_back({ S, cur });
        }
      }
    }
    
    void harvestRightPath(Production* start, const Production* stop, vector<pair<string, Production*>>& partition) {
      NonTerminal* cur = start->parent();
      Production* prev = start;
      while (cur != stop) {
        string S = "";
        while (S.size() < x_ && cur != stop) {
          if (cur->left() != prev) {
            assert(!cur->left()->associatedString.empty());
            assert(cur->left()->associatedString.length() < x_);
            
            assert(false);
            
            stringstream ss;
            ss << cur->left()->associatedString << S;
            S = ss.str();
          }
          
          prev = cur;
          cur = cur->parent();
        }
        
        if (!S.empty()) {
          partition.push_back({ S, cur });
        }
      }
    }
    
    void harvestIntermediateNodes(Production* a, Production* b, vector<pair<string, Production*>>& partition) {
      const NonTerminal* stop = lca(a, b);
      
      harvestLeftPath(a, stop, partition);
      harvestRightPath(b, stop, partition);
    }
    
    static vector<pair<string, Production*>> constructPartition(const SLP& slp, uint64_t x) {
      Partitioner partitioner(x);
      
      // Find key vertices
      partitioner.keyVertices.reserve(slp.derivedLength() / x);
      slp.root()->accept(&partitioner);
      
      // Find vertices in between
      vector<pair<string, Production*>> partition;
      partition.reserve(slp.derivedLength() / x);
      partitioner.harvestRightPath(partitioner.keyVertices[0], slp.root(), partition);
      for (int i = 0; i < partitioner.keyVertices.size() - 1; i++) {
        partition.push_back({ partitioner.keyVertices[i]->associatedString, partitioner.keyVertices[i] });
        partitioner.harvestIntermediateNodes(partitioner.keyVertices[i], partitioner.keyVertices[i + 1], partition);
      }
      auto lastKeyVertex = partitioner.keyVertices[partitioner.keyVertices.size() - 1];
      partition.push_back({ lastKeyVertex->associatedString, lastKeyVertex });
      partitioner.harvestLeftPath(lastKeyVertex, slp.root(), partition);
      
      return partition;
    }
    
  private:
    Partitioner(uint64_t x) : x_(x) { }
    
    const uint64_t x_;
    vector<Production*> keyVertices;
  };
}
