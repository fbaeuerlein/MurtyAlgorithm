/*
 * AssociationFinderMiller.h
 *
 *  Created on: 15.08.2013
 *      Author: fb
 */

#ifndef MILLER_H_
#define MILLER_H_

#include "AuctionAlgorithm.h"
#include <queue>


template<typename Scalar = double>
class AssociationFinderMiller {
public:

    typedef Eigen::Matrix<Scalar, -1, -1> WeightMatrix;
    typedef Eigen::Matrix<size_t, -1, -1> AssignmentMatrix;

    /**
     * a partition represents an assignment matrix with it's
     * weight matrix
     * see Murty's algorithm for details
     */
    class Partition
    {
    public:
        Partition() : value(0)
        {
            w = WeightMatrix::Zero(w.rows(), w.cols());
        }

        Partition(const Edges & edges, const WeightMatrix & w, const Scalar v) :
        edges(edges), w(w), value(v)
        {}

        Edges edges;
        WeightMatrix w;
        Scalar value;
    };

    struct ComparePartition: std::binary_function<Partition,Partition,bool>
    {
      bool operator()(const Partition & lhs, const Partition & rhs) const
      {
          return ( lhs.value < rhs.value );
      }
    };

    /**
     * list of partitions
     */
    typedef typename std::vector<Partition> Partitions;

    /**
     * sum up values of edges, i.e. objective function value
     * @param edges
     * @return
     */
    static Scalar objectiveFunctionValue(const Edges & edges, const size_t toX, const size_t toY )
    {
        Scalar v = 0;
        for ( auto & e : edges )
        {
            if ( e.y < toY)
                v += e.v;
        }
        return v;
    }

    static void edgesToAssignmentMatrix(AssignmentMatrix & m, const Edges & edges)
    {
        for(auto & e : edges )
            m(e.x, e.y) = 1;
    }

    static typename std::vector<Edges>  buildMBestAssignmentsAuction(const WeightMatrix & w, const size_t toX, const size_t toY, const size_t mBest = 5, const APSolvingMode mode = MAXIMIZATION)
    {
        const size_t rows = w.rows(), cols = w.cols();

        assert( rows != 0 && cols != 0 && cols >= rows );

        WeightMatrix m = w;

        typename std::vector<Edges> resultingEdges;

        // special case if rows = cols = 1
        if ( cols == 1 && rows == 1 )
        {
            if (m(0, 0) == 0 ) return resultingEdges;

            Edges edges;
            edges.push_back(Edge<Scalar>(0, 0, m(0, 0)));
            resultingEdges.push_back(edges);
            return resultingEdges;
        }

        size_t kBest = 0;

        const size_t maxComb = ( m.rows() > m.cols() ) ? m.rows() : m.cols();
        // if rows! < mBest ...
        switch (maxComb)
        {
            case 1 : kBest = 1; break;
            case 2 : kBest = 2; break;
            case 3 : kBest = 6; break;
            case 4 : kBest = 24; break;
            default: kBest = mBest; break;
        }
        if ( mBest < kBest ) kBest = mBest;

#ifdef __ASSOCIATON_FINDER_DEBUG
        std::cout << "try to solve: \n" << m << std::endl;
#endif
        Edges edges = Auction<Scalar>::solve(m); // make initial (best) assignment

        // sort edges by row
        std::sort(edges.begin(), edges.end(), [](const Edge<Scalar> & e1, const Edge<Scalar> & e2) {return e1.x < e2.x;});

        // initial partition, i.e. best solution
        Partition init(edges, m, objectiveFunctionValue(edges, toX, toY));

#ifdef __ASSOCIATON_FINDER_DEBUG
        {
            AssignmentMatrix a = AssignmentMatrix::Zero(m.rows(), m.cols());
            edgesToAssignmentMatrix(a, edges);

            std::cout << "generated initial assignment with value " << init.value
            << " and matrix =\n" << a << std::endl;
        }
#endif

        typedef std::priority_queue<Partition, std::vector<Partition>, ComparePartition > PartitionsPriorityQueue;

        // create answer-list with initial partition
        PartitionsPriorityQueue priorityQueue, answerList;
        priorityQueue.push(init);

        // assume values between 0 and 1 !
        const Scalar lockingValue = 0.;

        while ( !priorityQueue.empty() && answerList.size() < kBest )
        {
            // take first element from queue
            Partition currentPartition = priorityQueue.top();
            priorityQueue.pop();

            answerList.push(currentPartition);

            // for all triplets in this solution
            for (size_t e = 0; e < currentPartition.edges.size(); ++e)
            {
                auto & triplet = currentPartition.edges[e];

                // omit invalid triplets
                if ( triplet.y >= toY ) continue;

                WeightMatrix P_ = currentPartition.w; // P' = P

                // exclude edge by setting weight in matrix to lockingValue -> NOT (x, y)
                P_(triplet.x, triplet.y) = lockingValue;

#ifdef __ASSOCIATON_FINDER_DEBUG
                for (auto & t : currentPartition.edges)
                {
                    if ( t.x == triplet.x && t.y == triplet.y )
                    std::cout << "NOT ";
                    std::cout << "(" << t.x << ", " << t.y << ") ";
                }
                std::cout << std::endl;
                std::cout << P_ << std::endl;
#endif

                // determine solution for changed matrix and create partition
                Edges S_ = Auction<Scalar>::solve(P_);

                if (S_.size() == P_.rows())// solution found?
                {
                    // sort edges by row
                    std::sort(S_.begin(), S_.end(), [](const Edge<Scalar> & e1, const Edge<Scalar> & e2) {return e1.x < e2.x;});

                    Partition newPartition(S_, P_, objectiveFunctionValue(S_, toX, toY));

                    // if S exists
                    priorityQueue.push(newPartition);// push back unpartitioned new partition
                }
                // remove all vertices that include row and column of current node
                // i.e. force using this edge
                for (size_t r = 0; r < currentPartition.w.rows(); ++r )
                    currentPartition.w(r, triplet.y) = lockingValue;

                for (size_t c = 0; c < currentPartition.w.cols(); ++c )
                    currentPartition.w(triplet.x, c) = lockingValue;

                // set edge back to original value
                currentPartition.w(triplet.x, triplet.y) = triplet.v = m(triplet.x, triplet.y);

#ifdef __ASSOCIATON_FINDER_DEBUG
                std::cout << "setting (" << triplet.x << ", " << triplet.y << ") back to " << triplet.v << std::endl;
#endif
            }

        }

        // create return list
        while( !answerList.empty() )
        {
            resultingEdges.push_back(answerList.top().edges);
            answerList.pop();
        }

        return resultingEdges;
    }

};}
#endif
