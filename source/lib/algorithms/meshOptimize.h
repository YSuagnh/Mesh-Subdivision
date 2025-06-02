#ifndef MESH_OPTIMIZE_H
#define MESH_OPTIMIZE_H

#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <typeinfo>
#include <list>
#include <queue>
#include <random>
#include <mutex>
#include "MeshTopology.h"
#include "solverBase.h"

namespace std
{
    template<>
    struct hash<meshTopology::Point>
    {
        typedef size_t result_type;
        typedef meshTopology::Point argument_type;
        result_type operator() (const argument_type &p) const
        {
            return hash<meshTopology::real_t>()(p.p[0]) ^
                hash<meshTopology::real_t>() (p.p[1]) ^
                hash<meshTopology::real_t>() (p.p[2]);
        }
    };
}

namespace MeshOptimize
{
using namespace meshTopology;

struct Coord
{
    SDFace *f;
    std::array<double, 3> ks;
    double dist;
};
std::ostream &operator<< (std::ostream &os, const Coord &c)
{
    std::cout << "Coord{ " << *c.f;
    std::cout << "coords{";
    std::cout << c.ks[0] << "," << c.ks[1] << "," << c.ks[2] << "}";
    std::cout << "dist:" << c.dist << "}";
    return os;
}

// typedef std::vector<std::pair<int, real_t>> baryCoordType;
typedef std::vector<std::pair<SDVertex*, real_t>> CoordType;
typedef std::vector<std::vector<std::pair<size_t, double>>> SpM;

auto Energy(const SDShape&, const SDAttrib&, std::unordered_set<Point>)-> double;
auto ProjectPoints(const std::unordered_set<Point>, const std::vector<SDFace*>&, std::unordered_map<Point, Coord>&)-> double;
auto GenerateLegalMove(SDEdge, SDShape&, SDAttrib&, std::vector<SDEdge>&)-> bool;
auto ImproveVertexPosition(const SDShape, const std::unordered_map<Point, Coord>, const SDAttrib, SDVertex*)-> void;
auto ImproveVertexPositionV2(const SDShape, const std::unordered_map<Point, Coord>, const SDAttrib, SDVertex*)-> void;
static auto meshOptimize(tinyobj::attrib_t&, std::vector<tinyobj::shape_t>&)-> void;
auto star(SDVertex*, std::set<SDEdge>&, SDShape&, SDAttrib&)-> void;
auto star(SDEdge, SDShape&, SDAttrib&)-> void;

float SPRING_CONST = 0.001;
float crep = 0.001;
float THRESHOLD = 0.01;
float FACE_ANGLE = 0.5;
std::vector<Point> X;   // given data points
std::unordered_map<SDFace*, std::unordered_set<Point>> fs2x;   // data points that are projected to a given face


class meshOptimSolver: public SolverBase
{
    const char *name = "Mesh Optimize";
public:
    meshOptimSolver() {
        addFltParam("SPRING_CONST", SPRING_CONST);
        addFltParam("crep", crep);
        addFltParam("improve threshold", THRESHOLD);
        addFltParam("face angle", FACE_ANGLE);
    }
    const char *getName() {return name;}
    void run(tinyobj::attrib_t &a, std::vector<tinyobj::shape_t> &s) {
        SPRING_CONST = getFltParam(0);
        crep = getFltParam(1); 
        THRESHOLD = getFltParam(2);
        FACE_ANGLE = getFltParam(3);

        meshOptimize(a, s);
    }


    // ref: "Mesh Optimizing", Hugues Hoppe et. al
    void meshOptimize(tinyobj::attrib_t &attrib, std::vector<tinyobj::shape_t> &shapes)
    {
        std::vector<SDShape> Ks;
        SDAttrib V0;

        createSDShapes(attrib, shapes, V0, Ks);
        retainNormalVector(V0, Ks);

        // set V0 as given data points X
        for (SDVertex *v: V0.vs)
        {
            Point x = Point(v->p);
            X.push_back(x);
        }


        int shapeNum = 0;
        for (int i = 0; i < Ks.size(); ++i)
        {
            SDShape &K = Ks[i];

#ifndef NDEBUG
        std::cout << K << std::endl;
#endif

            std::cout << "shape" << shapeNum << ":\n";
            for (SDFace *fs: K.fs)
                for (SDVertex *v: fs->verts)
                    fs2x[fs].insert(v->p);

            std::vector<SDEdge> orders(K.es.begin(), K.es.end());
            std::random_device rd;
            std::mt19937 gr(rd());
            std::shuffle(orders.begin(), orders.end(), gr);

            // solve the outer minimization problem
            int edgeNum = 0;
            std::cout << "processed/need to be processed" << ":\n";
            while (edgeNum < orders.size())
            {
                auto eIt = K.es.find(orders[edgeNum++]);
                while (eIt == K.es.end())
                    eIt = K.es.find(orders[edgeNum++]);
                SDEdge edge = *K.es.find(*eIt);
                GenerateLegalMove(edge, K, V0, orders);

                if (edgeNum % 200 == 0)
                {
                    std::cout << edgeNum << "/" << orders.size() << "\n";

                    updateMtx.lock();
                    update = true;
                    sourceMtx.lock();
                    for (std::vector<SDVertex*>::size_type i = 0; i < V0.vs.size(); ++i)
                        V0.vs[i]->index = i;
                    retainNormalVector(V0, Ks);
                    createObjShapes(V0, Ks, attrib, shapes);
                    sourceMtx.unlock();
                    updateMtx.unlock();
                }
            }
            ++shapeNum;
        }

        finMtx.lock();
        hasFinished = true;
        sourceMtx.lock();
        for (std::vector<SDVertex*>::size_type i = 0; i < V0.vs.size(); ++i)
            V0.vs[i]->index = i;
        retainNormalVector(V0, Ks);
        createObjShapes(V0, Ks, attrib, shapes);
        sourceMtx.unlock();
        finMtx.unlock();
    }
};

double Energy(const SDShape &shape, const SDAttrib &V, std::unordered_set<Point> xs)
{
    double E = 0;
    std::unordered_map<Point, Coord> a;

#ifndef NDEBUG
    std::cout << "===Compute energy===\n";
    std::cout << shape << std::endl;
    std::cout << V << std::endl;
    for (Point x: xs)
        std::cout << x;
    std::cout << "\n";
#endif

    E += ProjectPoints(xs, shape.fs, a);    // E dist
    E += crep * V.vs.size();    // E rep
    for (const SDEdge e: shape.es)
    { // E spring
        Point diff = e[0]->p - e[1]->p;
        E += SPRING_CONST * diff * diff;
    }

#ifndef NDEBUG
    std::cout << "energy: " << E << "\n";
    std::cout << "======\n";
#endif

    return E;
}

void OptimizeVertexPosition(SDVertex *v, const SDShape &shape, 
    const SDAttrib &attrib, const std::unordered_set<Point> &X)
{
#ifndef NDEBUG
    std::cout << "before optimize: ";
    std::cout << shape << "\n";
    std::cout << "control points:\n";
    for (Point p: X)
        std::cout << p << "\n";
#endif

    // use the maximum point method
    for (int i = 0; i <18; ++i)
    {
        std::unordered_map<Point, Coord> B;
        ProjectPoints(X, shape.fs, B);
        ImproveVertexPosition(shape, B, attrib, v);
    }

    // use the conjugate gradient method
/*
    for (int i = 0; i < 12; ++i)
    {
        std::unordered_map<Point, Coord> B;
        ProjectPoints(X, shape.fs, B);
        ImproveVertexPositionV2(shape, B, attrib, v);
    }
    */

#ifndef NDEBUG
    std::cout << "after optimize: ";
    std::cout << shape << "\n";
#endif
}

/* ref: Philip J. Schneider and David H. Eberly, "geometric tools for computer graphics", chapter 10.3.2
    region notation:
             2
            \ |s
             \|
              |\
             3| \
              |  \
              | 0 \ 1
            __|____\_____t
            4 |  5  \ 6
*/
double ProjectPoints(const std::unordered_set<Point> X, const std::vector<SDFace*> &fs, 
    std::unordered_map<Point, Coord> &coords)
{
    static const int modulo[5]{0, 1, 2, 0, 1};

    double totalDist = 0;
    for (const auto &x: X)
    {
        Coord minCoord{nullptr, 0, 0, 0, MAXFLOAT};

        // for every SDFace f, calculate dist between x and f
        for (SDFace *face: fs)
        {
            Point e0 = face->verts[1]->p - face->verts[0]->p;
            Point e1 = face->verts[2]->p - face->verts[0]->p;
            Point BP = face->verts[0]->p - x;
            double a = e0 * e0, b = e0 * e1, c = e1 * e1;
            double d = BP * e0, e = BP * e1, f = BP * BP;

            auto dist = [=](double s, double t)->double
            {
                Point a = BP + s*e0 + t*e1;
                return a * a;
            };

            // get global minimum
            double denoG = a*c - b*b;
            double sbar = (b*e-c*d)/denoG;
            double tbar = (b*d-a*e)/denoG;

            // examine the region
            if (sbar >= 0)
            {
                if (tbar >= 0)
                {
                    if (sbar + tbar <= 1)
                    { // region 0
                        double alpha = 1 - sbar - tbar;
                        if (dist(sbar, tbar) < minCoord.dist)
                            minCoord = Coord{face, alpha, sbar, tbar, dist(sbar, tbar)};
                    } else
                    { // region 1
                        double deno1 = a - 2*b +c;
                        double t = (a-b+d-e)/deno1;
                        if (t < 0)
                            t = 0;
                        else if(t >1)
                            t = 1;
                        if (dist(1-t, t) < minCoord.dist)
                            minCoord = Coord{face, 0, 1-t, t, dist(1-t, t)};
                    }
                } else
                {
                    if (sbar + tbar < 1)
                    { // region 3
                        double s = -d / a;
                        if (s < 0)
                            s = 0;
                        else if (s > 1)
                            s = 1;
                        if (dist(s, 0) < minCoord.dist)
                            minCoord = Coord{face, 1-s, s, 0, dist(s, 0)};
                    } else
                    { // region 2
                        double mfs = -(a + d), ft = b + e;
                        double s = 0, t = 0;
                        if (mfs < 0)
                        { // t=0
                            t = 0, s = -d/a;
                            if (s < 0)
                                s = 0;
                            else if (s > 1)
                                assert(0);
                        } else if (mfs + ft < 0)
                        { // s+t=1
                            t = (a-b+d-e)/(a-2*b+c);
                            if (t>1)
                                t = 1;
                            else if (t < 0)
                                assert(0);
                        } else 
                        {
                            s = 1;
                        }
                        if (dist(s, t) < minCoord.dist)
                            minCoord = Coord{face, 1-s-t, s, t, dist(s, t)};
                    }
                }
            } else
            {
                if (tbar < 0)
                { //region 4
                    double fs = d, ft = e;
                    if (fs < 0)
                    { // t=0
                        double s = -d/a;
                        if (s > 1)
                            s = 1;
                        else if (s < 0)
                            assert(0);
                        if (dist(s, 0) < minCoord.dist)
                            minCoord = Coord{face, 1-s, s, 0, dist(s, 0)};
                    } else if (ft < 0)
                    {
                        double t = -e/c;
                        if (t > 1)
                            t = 1;
                        else if (t < 0)
                            assert(0);
                        if (dist(0, t) < minCoord.dist)
                            minCoord = Coord{face, 1-t, 0, t, dist(0, t)};
                    } else
                    {
                        if (dist(0, 0) < minCoord.dist)
                            minCoord = Coord{face, 1, 0, 0, dist(0, 0)};
                    }
                } else if (sbar + tbar < 1)
                {// region 5
                    double t = -e/c;
                    if (t < 0)
                        t = 0;
                    else if (t > 1)
                        t = 1;
                    if (dist(0, t) < minCoord.dist)
                        minCoord = Coord{face, 1-t, 0, t, dist(0, t)};
                } else
                { //region 6
                    double mft = -(c+e), fs = b+d;
                    if (mft < 0)
                    { // s=0
                        double t = -e / c;
                        if (t < 0)
                            t = 0;
                        else if (t > 1)
                            assert(0);
                        if (dist(0, t) < minCoord.dist)
                            minCoord = Coord{face, 1-t, 0, t, dist(0, t)};
                    } else if (fs + mft < 0)
                    { // s+t=1
                        double t = (a-b+d-e)/(a-2*b+c);
                        if (t < 0)
                            t = 0;
                        else if (t > 1)
                            t = 1;
                        if (dist(1-t, t) < minCoord.dist)
                            minCoord = Coord{face, 0, 1-t, t, dist(1-t, t)};
                    }
                }
            }
        }
        if (minCoord.f)
        {
            totalDist += minCoord.dist;
            coords.emplace(x, minCoord);
            fs2x[minCoord.f].insert(x);
        }
    }
    return totalDist;
}


// another way to do the same work with "conjugate gradient method"
void ImproveVertexPosition(const SDShape shape, const std::unordered_map<Point, Coord> B, 
    const SDAttrib attrib, SDVertex *v)
{
    real_t k = 0;
    Point sumA, sumBary;
    std::vector<SDVertex*> vs;
    for (SDEdge e: shape.es)
    {
        if (e[0] != v)
            vs.push_back(e[0]);
        else
            vs.push_back(e[1]);
    }
    for (SDVertex *ov: vs)
        sumA = sumA + ov->p;
    for (const auto& [x, coord]: B)
    {
        Point cur;
        real_t gamma;
        for (int i: {0, 1, 2})
        {
            if (coord.f->verts[i] != v)
                cur = cur + coord.ks[i]* coord.f->verts[i]->p;
            else
                gamma = coord.ks[i];
        }
        cur = cur - x;
        cur = cur * gamma;
        sumBary = sumBary + cur;
        k += gamma * gamma;
    }
    k += shape.es.size() * SPRING_CONST;
    v->p = (SPRING_CONST * sumA - sumBary) / k;
}

/* input: A, x
* return: Ax
* sparse matrix-vector multiplication(SpMV).
* y = Ax, where matrix A is sparse and vector x is dense.
*/
inline std::vector<double> SpMV(const SpM &A, std::vector<double> &x)
{
    std::vector<double> y;
    for (const auto &row: A)
    {
        double res = 0;
        for (const auto &[i, val]: row)
            res += val * x[i];
        y.push_back(res);
    }
    return y;
}

/* input: A
* return: A^T
*/
inline SpM SpMInverse(SpM A, size_t cols)
{
    SpM AT;
    int curCols = 0;
    for (int i = 0; auto &row: A)
    {
        for (auto &[j, val]: row)
        {
            while (AT.size() <= j)
                AT.emplace_back();
            AT[j].emplace_back(i, val);
        }
        ++i;
    }
    while (AT.size() < cols)
        AT.emplace_back();
    return AT;
}

// vector add
// c = a + scale * b, element-wise, size of a and b must be the same.
inline std::vector<double> VecAdd(const std::vector<double> &a, 
    const std::vector<double> &b, double scale=1)
{
    if (a.size() != b.size()) std::cout << "size wrong\n";
    std::vector<double> c;
    for (std::vector<double>::size_type i = 0; i < a.size(); ++i)
        c.push_back(a[i] + scale * b[i]);
    return c;
}

// vector sub
// c = a - scale*b, element-wise, size of a and b must be the same.
inline std::vector<double> VecSub(const std::vector<double> &a, 
    const std::vector<double> &b, const double scale=1)
{
    if (a.size() != b.size()) std::cout << "size wrong\n";
    std::vector<double> c;
    for (std::vector<double>::size_type i = 0; i < a.size(); ++i)
        c.push_back(a[i] - scale * b[i]);
    return c;
}

// inner-product of two vectors
// return a*b, where * is inner-product
inline double VecDot(const std::vector<double> &a, const std::vector<double> &b)
{
    if (a.size() != b.size()) std::cout << "size wrong\n";
    double res = 0;
    for (std::vector<double>::size_type i = 0; i < a.size(); ++i)
        res += a[i] * b[i];
    return res;
}

// ref: chapter 11.3.9(CGNR) "Matrix Computations(4th edition)", Gene H. Golub, Charles F. Van Loan
void ImproveVertexPositionV2(const SDShape K, const std::unordered_map<Point, Coord> B,
    const SDAttrib attrib, SDVertex *v)
{
    // find all related vertices
    std::set<SDVertex*> vSet;
    std::vector<SDVertex*> vs;
    std::vector<Point> ps;
    for (SDFace *f: K.fs)
        for (SDVertex *v: f->verts)
            vSet.insert(v);
    vs.assign(vSet.begin(), vSet.end());

    // indexing the vertex
    std::map<SDVertex*, size_t> vIndex;
    for (size_t i = 0; auto *v: vs)
        vIndex[v] = i++;

    // construct sparse matrix A and vector d
    SpM A;
    for (auto &[p, coord]: B)
    {
        ps.push_back(p);
        A.emplace_back();
        auto &row = A.back();
        for (int i = 0; double k: coord.ks)
            row.emplace_back(vIndex[coord.f->verts[i++]], k);

    }
    for (SDEdge e: K.es)
    {
        A.emplace_back();
        auto &row = A.back();
        row.emplace_back(vIndex[e[0]], sqrt(SPRING_CONST));
        row.emplace_back(vIndex[e[1]], -sqrt(SPRING_CONST));
    }

    Point newP;
    for (int i: {0, 1, 2})
    {
        // build vector d
        std::vector<double> d;
        for (auto &p: ps)
            d.push_back(p[i]);
        for (auto e: K.es)
            d.push_back(0);
        // build vector x
        std::vector<double> x;
        for (SDVertex *v: vs)
            x.push_back(v->p[i]);
        
        // begin the main algorithm
        SpM AT = SpMInverse(A, vs.size());
        auto r = VecSub(d, SpMV(A, x));
        auto z = SpMV(AT, r);
        auto p = z;
        int ii = 0;
        while (VecDot(z, z) > 1e-3)
        {
            double miu = VecDot(z, z) / (VecDot(SpMV(A, p), SpMV(A, p)));
            x = VecAdd(x, p, miu);
            r = VecSub(r, SpMV(A, p), miu);
            double zcdot = VecDot(z, z);
            z = SpMV(AT, r);
            double tau = VecDot(z, z) / zcdot;
            p = VecAdd(z, p, tau);
        }

        newP[i] = x[vIndex[v]];
    }
    v->p = newP;
}

bool checkDiheralAngle(const std::set<SDEdge> &es)
{
#ifndef NDEBUG
    std::cout << "===check diheral angle===\n";
#endif

    for (SDEdge e: es)
    {
        if (!e.fs[1]) continue;
        // check diheral angle between neighbor faces
        // their normal vectors are unit

#ifndef NDEBUG
        std::cout << *e.fs[0] << e.fs[0]->vns[0]->p << " and ";
        std::cout << *e.fs[1] << e.fs[1]->vns[0]->p << "\n";
#endif

        if (e.fs[0]->vns[0]->p * e.fs[1]->vns[0]->p < FACE_ANGLE)
            return false;
    }
    return true;
}

// recover neighbors, edges, startFace info in SDShape.
void localRecover(SDShape &shape)
{
    int modulo[5]{0, 1, 2, 0, 1};
    std::set<SDEdge> es;
    std::set<SDEdge> resEs;
    for (SDFace *f: shape.fs)
    {
        for (int i: {0, 1, 2})
        {
            f->verts[i]->startFace = f;
            SDEdge e(f->verts[modulo[i]], f->verts[modulo[i+1]]);
            if (shape.es.find(e) == shape.es.end())
                continue;
            auto place = es.find(e);
            if (place == es.end())
            {
                e.fs[0] = f;
                e.f0edgeNum = i;
                es.insert(e);
            } else
            {
                SDEdge prev = *place;
                prev.fs[1] = f;
                prev.f1edgeNum = i;
                resEs.insert(prev);
                
                prev.fs[0]->neighborFs[prev.f0edgeNum] = f;
                f->neighborFs[i] = prev.fs[0];
            }
        }
    }
    for (SDEdge e: es)
        resEs.insert(e);
    shape.es.swap(resEs);

#ifndef NDEBUG
    std::cout << "*after recover*\n";
    std::cout << shape << "\n";
    for (SDFace *f: shape.fs)
    {
        std::cout << "neighbor of " << *f << ":\n";
        for (SDFace *nf: f->neighborFs)
        {
            if (nf)
                std::cout << *nf;
        }
        std::cout << std::endl;
    }
    for (SDEdge e: shape.es)
    {
        std::cout << "neighbor face of " << e << ":\n";
        std::cout << *e.fs[0];
        if (e.fs[1])
            std::cout << *e.fs[1];
        std::cout << std::endl;
    }
#endif
}

void star(SDVertex *v, std::set<SDEdge> &es, SDShape &outSdshape, SDAttrib &outAttrib)
{
    std::vector<SDVertex*> vs;
    findFsByVtx(v, std::back_inserter(outSdshape.fs));
    findVtxNeighbors(v, std::back_inserter(vs));
    for (SDVertex *vtx: vs)
    {
        outSdshape.es.insert(*es.find(SDEdge(vtx, v)));
    }
    outAttrib.vs.push_back(v);
}

void star(SDEdge e, SDShape &sdshape, SDAttrib &attrib)
{
    sdshape.fs.push_back(e.fs[0]);
    sdshape.fs.push_back(e.fs[1]);
    sdshape.es.insert(e);
}

bool GenerateLegalMove(SDEdge e, SDShape &K, SDAttrib &V, std::vector<SDEdge> &orders)
{
    bool cond = false;
    static const int modulo[5]{0, 1, 2, 0, 1};

#ifndef NDEBUG
    std::cout << "///GenerateLegalMove///\n";
    std::cout << "current shape: \n";
    std::cout << K << "\n";
    std::cout << "selected edge: ";
    std::cout << e << std::endl;
#endif

    // SITUATION 1: collapsing edge
    bool cond0=true, cond1=true, cond2=false;
    std::set<SDVertex*> vs0, vs1;
    findVtxNeighbors(e[0], std::inserter(vs0, vs0.begin()));
    findVtxNeighbors(e[1], std::inserter(vs1, vs1.begin()));
    for (SDVertex *v: vs0)
    {
        /*
        cond0 must be checked. Consider this structure:
                             ______________ 
                            |\*            |              
                            | \  *   face3 |
                            |  \    *      |           
                            |   \ face0*   |
                         e0 |face\________*.
                            |    /        *|
                            |   / face1*   |
                            |  /    *      |
                            | /  *   face2 | 
                            |/*____________|
        If we collapse e0, a hanging face will be created and
        program will crash or stuck.
        */
        if (v == e[0] || v == e[1]) continue;
        if (vs1.find(v) != vs1.end() && v != e.fs[0]->otherVertex(&e))
        {
            if (e.fs[1] && v != e.fs[1]->otherVertex(&e))
                cond0 = false;
            else if (!e.fs[1])
                cond0 = false;
        }
    }
    if (isBoundary(e[0]) && isBoundary(e[1]))
    {
        if (!e.boundary())
            cond1 = false;
    }
    if (isBoundary(e[0]) || isBoundary(e[1]))
    {
        /*
        This condition ensures that program will not try to do edge collapse
        on a tetrahedron or triangle: these will cause dimensionality reduction 
        problem for model.
        */
        std::set<SDVertex*> vs;
        for (SDFace *f: K.fs)
        {
            for (int i: {0, 1, 2})
                vs.insert(f->verts[i]);
            if (vs.size() > 3)
            {
                cond2 = true;
                break;
            }
        }
    } else
    {
        std::set<SDVertex*> vs;
        for (SDFace *f: K.fs)
        {
            for (int i: {0, 1, 2})
                vs.insert(f->verts[i]);
            if (vs.size() > 4)
            {
                cond2 = true;
                break;
            }
        }
    }
    cond = cond0 && cond1 && cond2;
    if (cond)
    {
#ifndef NDEBUG
        std::cout << "===Move: edge collapse===\n";
#endif
        // select loacl simplicial complex(star({i}, K) U star({j}, K))
        SDShape oldShape;
        SDAttrib oldAttrib;
        std::set<SDFace*> fs;
        std::set<SDEdge> es;
        fs.insert(e.fs[0]);
        findFsByVtx(e[0], std::inserter(fs, fs.begin()));
        findFsByVtx(e[1], std::inserter(fs, fs.begin()));
        for (SDVertex *v: vs0)
            es.emplace(e[0], v);
        for (SDVertex *v: vs1) 
        {
            es.emplace(e[1], v);
            vs0.insert(v);
        }
            
        oldShape.fs.assign(fs.begin(), fs.end());
        oldShape.es.swap(es);
        oldAttrib.vs.assign(vs0.begin(), vs0.end());

#ifndef NDEBUG
        std::cout << "energy before imrove:\n";
#endif
        // compute energy on selected simplicial complex
        std::unordered_set<Point> xs;
        for (SDFace *f: oldShape.fs)
            for (Point p: fs2x[f])
                xs.insert(p);
        double oldEnergy = Energy(oldShape, oldAttrib, xs);

        std::vector<SDShape> newShapes;
        std::vector<SDAttrib> newAttribs;
        std::vector<SDFace*> fsSurroundA, fsSurroundB;
        findFsByVtx(e[0], std::back_inserter(fsSurroundA));
        findFsByVtx(e[1], std::back_inserter(fsSurroundB));
        // possible positions
        Point initPos[3]{(e[0]->p+e[1]->p)/2, e[0]->p, e[1]->p};
        SDVertex *newVs[3]{0, 0, 0};
        int best = 0;
        float bestEnergy = MAXFLOAT;
        for (int i: {0, 1, 2})
        {
            newShapes.emplace_back();
            newAttribs.emplace_back();
            // construct new simplicial complex
            SDShape &newShape = newShapes.back();
            SDAttrib &newAttrib = newAttribs.back();
            // possible initial value of new vertex
            newVs[i] = new SDVertex(initPos[i]);
            SDVertex *newV = newVs[i];
            // create new sdfaces
            auto createNewFace = [&](SDFace *fs, int i)-> void
            {
                if (fs != e.fs[0] && fs != e.fs[1])
                {
                    int j = fs->getVertIndex(e[i]);
                    for (SDFace *f: newShape.fs)
                        if ((f->verts[1]==(*fs)[modulo[j+1]] && f->verts[2]==(*fs)[modulo[j+2]]) ||
                            (f->verts[1]==(*fs)[modulo[j+2]] && f->verts[2]==(*fs)[modulo[j+1]]))
                                return;
                    SDFace *newF = new SDFace(newV, (*fs)[modulo[j+1]], (*fs)[modulo[j+2]]);
                    newShape.fs.push_back(newF);
                    // newF->neighborFs[1] = fs->neighborFs[modulo[j+1]];
                    newF->index = fs->index;
                }
            };
            std::for_each(fsSurroundA.begin(), fsSurroundA.end(), 
                std::bind(createNewFace, std::placeholders::_1, 0));
            std::for_each(fsSurroundB.begin(), fsSurroundB.end(), 
                std::bind(createNewFace, std::placeholders::_1, 1));
            retainNormalVector(newAttrib, newShape);
            // create new sdedge
            std::set<SDEdge> newEs;
            for (SDFace *f: newShape.fs)
            {
                for (int i:{0, 2})
                {
                    SDEdge e = SDEdge(f->verts[i], f->verts[modulo[i+1]]);
                    if (auto place = newEs.find(e); place == newEs.end())
                    {
                        e.fs[0] = f;
                        e.f0edgeNum = i;
                        newEs.insert(e);
                    } else
                    {
                        e.fs[0] = place->fs[0];
                        e.fs[1] = f;
                        e.f0edgeNum = place->f0edgeNum;
                        e.f1edgeNum = i;

                        e.fs[0]->neighborFs[e.f0edgeNum] = f;
                        f->neighborFs[i] = e.fs[0];
                        newShape.es.insert(e);
                        newEs.erase(place);
                    }
                }
            }
            for (SDEdge e: newEs)
                newShape.es.insert(e);
            // select sdvertex
            newAttrib.vs.push_back(newV);

            // improve vertex position
            std::unordered_set<Point> xcur;
            for (SDFace *f: oldShape.fs)
            {
                for (Point p: fs2x[f])
                    xcur.insert(p);
            }
            OptimizeVertexPosition(newV, newShape, newAttrib, xcur);

#ifndef NDEBUG
            std::cout << "energy after optimized\n" ;
            Energy(newShape, newAttrib, xs);
#endif
            
            if (float cur = Energy(newShape, newAttrib, xs); cur < bestEnergy)
            {
                bestEnergy = cur;
                best = i;
            }
        }
        retainNormalVector(newAttribs[best], newShapes[best]);

        SDShape reShape;
        std::set<SDFace*> fset(newShapes[best].fs.begin(), newShapes[best].fs.end());
        for (auto fs: oldShape.fs)
        {
            for (SDFace *nf: fs->neighborFs)
            {
                if (nf && std::find(oldShape.fs.begin(), oldShape.fs.end(), nf) == oldShape.fs.end())
                    fset.insert(nf);
            }
        }
        for (auto fs: newShapes[best].fs)
        {
            for (int i: {0, 1, 2})
                reShape.es.emplace(fs->verts[i], fs->verts[modulo[i+1]]);
        }
        reShape.fs.assign(fset.begin(), fset.end());
        localRecover(reShape);

        if (oldEnergy - bestEnergy >= oldEnergy * THRESHOLD &&
            checkDiheralAngle(reShape.es))
        {   
#ifndef NDEBUG
            std::cout << "accept move\n";
#endif

            for (SDFace *f: oldShape.fs)
                K.fs.erase(std::find(K.fs.begin(), K.fs.end(), f));
            for (SDFace *f: newShapes[best].fs)
                K.fs.push_back(f);
            for (SDVertex *v: e.verts)
                V.vs.erase(std::find(V.vs.begin(), V.vs.end(), v));
            V.vs.push_back(newVs[best]);

            SDAttrib erAttrib;
            std::set<SDVertex*> vs;
            for (SDFace *f: oldShape.fs)
                for (SDVertex *v: f->verts)
                    vs.insert(v);
            erAttrib.vs.assign(vs.begin(), vs.end());

            for (SDEdge e: oldShape.es)
                K.es.erase(e);
            for (SDEdge e: reShape.es)
            {
                K.es.erase(e);
                K.es.insert(e);
            }

            for (SDEdge newE: newShapes[best].es)
                orders.push_back(newE);

            // delete e[0];
            // delete e[1];
        } else
        {
#ifndef NDEBUG
            std::cout << "reject move\n";
#endif
            // rebuild neighbors and edges information
            SDShape reShape0;
            reShape0.es.swap(oldShape.es);
            std::set<SDFace*> fset0(oldShape.fs.begin(), oldShape.fs.end());
            for (int i = 0; auto fs: {fsSurroundA, fsSurroundB})
            {
                for (SDFace *f: fs)
                {
                    int index = f->getVertIndex(e[i]);
                    SDFace *nf = f->otherNeighbor(e[i]);
                    if (nf) fset0.insert(nf);
                    reShape0.es.emplace(f->verts[modulo[index+1]], f->verts[modulo[index+2]]);
                }
                ++i;
            }
            reShape0.fs.assign(fset0.begin(), fset0.end());
            localRecover(reShape0);
        
            /*
            for (SDFace *fp: newShapes[best].fs)
            {
                delete fp;
            }
            delete newVs[best];*/

        }

        /*
        for (int i = 0; i < 3; ++i)
        {
            if (i == best) continue;
            for (SDFace *fp: newShapes[i].fs)
            {
                delete fp;
            }
            delete newVs[i];
        }*/

        return false;
    }
    
    // SITUATION 2: edge swap
    SDVertex *v0 = e.fs[0]->otherVertex(&e), *v1 = nullptr;
    if (!e.boundary())
    {
        v1 = e.fs[1]->otherVertex(&e);
        SDEdge newE(v0, v1);
        if (K.es.find(newE) == K.es.end())
            cond = true;
    }
    if (cond)
    {
#ifndef NDEBUG
        std::cout << "===Move: edge swap===\n";
#endif

        Point oldPos[2]{v0->p, v1->p};
        double oldEnergies[2], energies[2];

        // construct new sdedge and sdface
        int otherIndex = e.fs[0]->getVertIndex(e.fs[0]->otherVertex(&e));
        SDVertex *orderVtx[2]{e.fs[0]->verts[modulo[otherIndex+1]], e.fs[0]->verts[modulo[otherIndex+2]]};
        SDVertex *otherVtx[2]{e.fs[0]->otherVertex(&e), e.fs[1]->otherVertex(&e)};
        SDEdge newEdge = SDEdge(otherVtx[0], otherVtx[1]);
        newEdge.fs[0] = new SDFace(otherVtx[0], orderVtx[0], otherVtx[1]);
        newEdge.fs[1] = new SDFace(otherVtx[0], otherVtx[1], orderVtx[1]);

        // consider two submesh
        std::vector<SDShape> subShapes;
        std::vector<SDAttrib> subAttribs;
        for (int i: {0, 1})
        {
            SDShape oldShape;
            SDAttrib oldAttrib;
            subShapes.emplace_back();
            subAttribs.emplace_back();
            SDShape &subShape = subShapes.back();
            SDAttrib &subAttrib = subAttribs.back();
            SDVertex *v = newEdge[i];

            star(v, K.es, oldShape, oldAttrib);

            std::unordered_set<Point> xs;
            for (SDFace *f: oldShape.fs)
                for (Point p: fs2x[f])
                    xs.insert(p);

            oldEnergies[i] = Energy(oldShape, oldAttrib, xs);

            // select meshs of v_h in K'
            subShape = oldShape;
            if (auto place = std::find(subShape.fs.begin(), subShape.fs.end(), e.fs[0]);
                place != subShape.fs.end())
                subShape.fs.erase(place);
            if (auto place = std::find(subShape.fs.begin(), subShape.fs.end(), e.fs[1]);
                place != subShape.fs.end())
                subShape.fs.erase(place);
            subShape.fs.push_back(newEdge.fs[0]);
            subShape.fs.push_back(newEdge.fs[1]);
            subShape.es.insert(newEdge);

            retainNormalVector(subAttrib, subShape);
            OptimizeVertexPosition(v, subShape, subAttrib, xs);
            energies[i] = Energy(subShape, subAttrib, xs);
        }

        // choose the best
        int best = energies[0] < energies[1]? 0: 1;
        double oldBestEnergy = oldEnergies[best];

        std::set<SDEdge> es = subShapes[0].es;
        es.insert(subShapes[1].es.begin(), subShapes[1].es.end());
        // check diheral angle
        if ((oldBestEnergy-energies[best]) > oldBestEnergy*THRESHOLD 
            && checkDiheralAngle(es))
        {
#ifndef NDEBUG
            std::cout << "accept move\n";
#endif

            SDShape reShape;
            std::set<SDFace*> reFs;
            std::set<SDEdge> reEs;
            for (int i: {0, 1})
            {
                reFs.insert(newEdge.fs[i]);
                for (SDFace *nf: e.fs[i]->neighborFs)
                    if (nf && nf != e.fs[1-i])
                        reFs.insert(nf);
                for (int j: {0, 1, 2})
                    reEs.emplace(newEdge.fs[i]->verts[j], newEdge.fs[i]->verts[modulo[j+1]]);
            }
            reShape.fs.assign(reFs.begin(), reFs.end());
            reShape.es.swap(reEs);
            localRecover(reShape);

            // update mesh
            for (SDFace *f: e.fs)
                K.fs.erase(std::find(K.fs.begin(), K.fs.end(), f));
            for (SDFace *f: newEdge.fs)
                K.fs.push_back(f);
            for (SDEdge e: reShape.es)
            {
                K.es.erase(e);
                K.es.insert(e);
            }
            // update orders
            K.es.erase(e);
            orders.push_back(newEdge);

            return true;
        } else
        {
#ifndef NDEBUG
            std::cout << "reject move\n";
#endif
            v0->p = oldPos[0];
            v1->p = oldPos[1];

            return false;
        }

    }

    // SITUATION 3: edge split
#ifndef NDEBUG
    std::cout << "===Move: edge split===\n";
#endif
    SDShape oldShape, newShape;
    SDAttrib oldAttrib, newAttrib;
    // select old simplicial complex
    if (e.fs[0]) oldShape.fs.push_back(e.fs[0]);
    if (e.fs[1]) oldShape.fs.push_back(e.fs[1]);
    oldShape.es.insert(e);
    
    std::unordered_set<Point> xs;
    for (SDFace *f: e.fs)
        for (Point p: fs2x[f])
            xs.insert(p);
    double oldEnergy = Energy(oldShape, oldAttrib, xs);

    // construct new simplicial complex
    // SDAttrib.vs
    Point newp((e[0]->p + e[1]->p)/2);
    SDVertex *newV = new SDVertex(newp);
    newAttrib.vs.push_back(newV);
    // SDShape.fs, SDShape.es
    for (int j: {0, 1})
    {
        SDFace *f = e.fs[j];
        if (!f) break;  // if e is boundary

        SDVertex *v = f->otherVertex(&e);
        int i = f->getVertIndex(v);
        newShape.fs.push_back(new SDFace(newV, v, (*f)[modulo[i+1]]));
        newShape.fs.push_back(new SDFace(newV, (*f)[modulo[i+2]], v));

        SDEdge e(newV, v);
        e.fs[0] = newShape.fs[2*j];
        e.fs[1] = newShape.fs[2*j+1];
        newShape.es.insert(e);
    }
    retainNormalVector(newAttrib, newShape);

    SDEdge e0(newV, (*newShape.fs[0])[2]);
    e0.fs[0] = newShape.fs[0];
    if(!e0.boundary()) e0.fs[1] = newShape.fs[3];
    newShape.es.insert(e0);

    SDEdge e1(newV, (*newShape.fs[1])[1]);
    e1.fs[0] = newShape.fs[1];
    if (!e.boundary()) e1.fs[1] = newShape.fs[2];
    newShape.es.insert(e1);

    // optimize vertex position
    std::unordered_set<Point> xcur;
    for (SDFace *f: oldShape.fs)
    {
        for (Point p: fs2x[f])
            xcur.insert(p);
    }
    OptimizeVertexPosition(newV, newShape, newAttrib, xcur);

    double newEnergy = Energy(newShape, newAttrib, xs);
    if (oldEnergy - newEnergy > oldEnergy * THRESHOLD &&
        checkDiheralAngle(newShape.es))
    {
#ifndef NDEBUG
        std::cout << "accept move\n";
#endif

        // rebuild mesh local information
        SDShape reShape = newShape;
        for (SDFace *f: oldShape.fs)
        {
            for (SDFace *nf: f->neighborFs)
            {
                if (nf && nf != e.fs[0] && nf != e.fs[1])
                    reShape.fs.push_back(nf);
            }
            for (int i: {0, 1, 2})
            {
                SDEdge curE = SDEdge(f->verts[i], f->verts[modulo[i+1]]);
                if (curE == e) 
                {
                    continue;
                }
                reShape.es.insert(curE);
            }
        }
        localRecover(reShape);
        for (SDEdge e: reShape.es)
        {
            K.es.erase(e);
            K.es.insert(e);

            if (newShape.es.find(e) != newShape.es.end())
                orders.push_back(e);
        }
        K.es.erase(e);

        for (SDFace *f: oldShape.fs)
            K.fs.erase(std::find(K.fs.begin(), K.fs.end(), f));
        for (SDFace *f: newShape.fs)
            K.fs.push_back(f);
        V.vs.push_back(newAttrib.vs[0]);

        return true;
    } else
    {
#ifndef NDEBUG
        std::cout << "reject move\n";
#endif
        return false;
    }
}

} // namespace: MeshOptimize

#endif