#ifndef COINCLIQUEEXTENDER_HPP
#define COINCLIQUEEXTENDER_HPP

class CoinCliqueList;
class CoinConflictGraph;

#include "CoinUtilsConfig.h"
#include <cstddef>

class COINUTILSLIB_EXPORT CoinCliqueExtender {
public:
    explicit CoinCliqueExtender(const CoinConflictGraph *cgraph, size_t extMethod = 4, const double *rc = NULL,
            const double maxRC = 100.0);
    ~CoinCliqueExtender();

    size_t extendClique(const size_t *clqIdxs, const size_t clqSize);

    size_t nCliques() const;
    const size_t* getClique(const size_t i) const;
    size_t getCliqueSize(const size_t i) const;

    void setMaxCandidates(const size_t maxCandidates);

private:
    void fillCandidates(const size_t *clqIdxs, const size_t clqSize);
    size_t randomExtension(const size_t *clqIdxs, const size_t clqSize);
    size_t greedySelection(const size_t *clqIdxs, const size_t clqSize, const double *costs);

    const CoinConflictGraph *cgraph_; //complete conflict graph
    CoinCliqueList *extendedCliques_; //stores the extended cliques

    /*Extending method:
     * 0 - no extension
     * 1 - random
     * 2 - max degree
     * 3 - max modified degree
     * 4 - reduced cost
     * 5 - reduced cost + modified degree
    */
    size_t extMethod_;

    //maximum size of the candidates list
    size_t maxCandidates_;

    //auxiliary arrays
    size_t *candidates_, nCandidates_; //used to temporarily store candidate vertices.
    size_t *newClique_, nNewClique_; //used to temporarily store a new extended clique.
    double *costs_;
    bool *iv_, *iv2_;

    //reduced cost parameters
    const double *rc_;
    double maxRC_;
};


#endif //COINCLIQUEEXTENDER_HPP
