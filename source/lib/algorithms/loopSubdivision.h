#ifndef LOOP_SUBDIVISION_H
#define LOOP_SUBDIVISION_H

#include <set>
#include <vector>
#include "MeshTopology.h"
#include "solverBase.h"

#ifndef NDEBUG
#include <iostream>
#endif


namespace LoopSubdivision 
{ // ref: https://pbr-book.org/3ed-2018/Shapes/Subdivision_Surfaces#SDVertex::regular
using namespace meshTopology;

static void loopSubdivision(tinyobj::attrib_t&, std::vector<tinyobj::shape_t>&);
inline real_t beta(int valence);
template <class V>
void computeWeightedPosition(const std::vector<V*> &comp, real_t beta, Point &resp);

class LoopSolver: public SolverBase
{
    const char *name = "Loop Subdivision";
public:
    const char *getName() {return name;}
    void run(tinyobj::attrib_t &attrib, std::vector<tinyobj::shape_t> &shapes) 
    {
        std::vector<SDShape> sdshapes;
        SDAttrib sdattrib;
        createSDShapes(attrib, shapes, sdattrib, sdshapes);

        loopSubdivision(sdattrib, sdshapes);

        finMtx.lock();
        hasFinished = true;
        sourceMtx.lock();
        retainNormalVector(sdattrib, sdshapes);
        createObjShapes(sdattrib, sdshapes, attrib, shapes);
        sourceMtx.unlock();
        finMtx.unlock();
    }

    /*
    * Perform Loop subdivision. Save the mesh after subdivision in sdattrib and sdshapes.
        normal vector and textcoord are not neeeded. In output, we just leave sdattrib.ns and sdattrib.ts empty.
    * input & output: sdattrib. sdattrib saves vertex, normal vector and textcoord of the mesh.
        sdshapes. sdshapes saves topology structure of the mesh.
    * ref: https://pbr-book.org/3ed-2018/Shapes/Subdivision_Surfaces
    */
    void loopSubdivision(SDAttrib &sdattrib, std::vector<SDShape> &sdshapes)
    {
        /* TODO: Complete this function*/
    }
};

inline real_t beta(int valence)
{
    if (valence == 3) return 3.f / 16.f;
    return 3.f / (8.f * valence);
}

template <class V>
void computeWeightedPosition(const std::vector<V*> &comp, real_t beta, Point &resp)
{
    for (int i = 0; i < 3; i++)
    {
        resp[i] *= (1 - comp.size() * beta);
        for (V *vertp: comp) 
            resp[i] += vertp->p[i] * beta;
    }
}


} // namespace: loopSubdivision

#endif