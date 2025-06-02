#ifndef ROOT3_SUBDIVISION_H
#define ROOT3_SUBDIVISION_H

#include <cmath>
#include "solverBase.h"
#include "MeshTopology.h"

#ifndef NDEBUG
#include <iostream>
#endif

namespace root3Subdivision
{
using namespace meshTopology;

static void root3Subdivision(tinyobj::attrib_t&, std::vector<tinyobj::shape_t>&);
template <class It>
inline void oneRingRule(const It begin, const It end, SDVertex *res, real_t w);
inline real_t alpha(int valence);

class Root3Solver: public SolverBase
{
    const char *name = "Root3 Subdivision";
public:
    const char *getName() {return name;}
    void run(tinyobj::attrib_t &attrib, std::vector<tinyobj::shape_t> &shapes) 
    {
        std::vector<SDShape> sdshapes;
        SDAttrib sdattrib;
        createSDShapes(attrib, shapes, sdattrib, sdshapes);

        root3Subdivision(sdattrib, sdshapes);

        finMtx.lock();
        hasFinished = true;
        sourceMtx.lock();
        retainNormalVector(sdattrib, sdshapes);
        createObjShapes(sdattrib, sdshapes, attrib, shapes);
        sourceMtx.unlock();
        finMtx.unlock();
    }

    /*
    * Perform Root3 subdivision. Save the mesh after subdivision in sdattrib and sdshapes.
        normal vector and textcoord are not needed. In output, we just leave sdattrib.ns and sdattrib.ts empty.
    * input & output: sdattrib. sdattrib saves vertex, normal vector and textcoord of the mesh.
        sdshapes. sdshapes saves topology structure of the mesh.
    * ref: https://dl.acm.org/doi/10.1145/344779.344835
    */
    void root3Subdivision(SDAttrib &sdattrib, std::vector<SDShape> &sdshapes)
    {
        /*TODO: complete this function*/
    }

};

template <class It>
inline void oneRingRule(const It begin, const It end, SDVertex *res, real_t w)
{
    real_t wres = 1 - res->valence * w;
    for (int i = 0; i < 3; i++)
        res->p[i] *= wres;
    for (It ele = begin; ele != end; ++ele)
    {
        for (int i: {0, 1, 2})
            res->p[i] += w * (*ele)->p[i];
    }
}

inline real_t alpha(int valence)
{
    switch (valence)
    {
    case 1:
        return 2.f / 9.f;
    case 2:
        return 1.f / 3.f;
    case 3:
        return 5.f / 27.f;
    case 4:
        return 1.f / 9.f;
    case 5:
        break;
    case 6:
        return 1.f / 18.f;    
    default:
        break;
    }
    const float pi = 3.141592653589;
    return (4.f-2.f*std::cos(2.f*pi/valence)) / (9.f*(float)valence);
}

}

#endif