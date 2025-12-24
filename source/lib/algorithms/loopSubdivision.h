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
        std::vector<SDVertex*> new_vs;

        // 1. Update Even Vertices
        for (SDVertex *v : sdattrib.vs) {
            Point p_new;
            if (isBoundary(v)) {
                std::vector<SDVertex*> neighbors;
                findVtxNeighbors(v, std::back_inserter(neighbors));
                if (neighbors.size() >= 2) {
                    Point p_neighbors = neighbors.front()->p + neighbors.back()->p;
                    p_new = 0.75 * v->p + 0.125 * p_neighbors;
                } else {
                    p_new = v->p; 
                }
            } else {
                std::vector<SDVertex*> neighbors;
                findVtxNeighbors(v, std::back_inserter(neighbors));
                int k = neighbors.size();
                real_t b = beta(k);
                Point sum_neighbors(0,0,0);
                for(auto* n : neighbors) sum_neighbors = sum_neighbors + n->p;
                p_new = (1.0 - k * b) * v->p + b * sum_neighbors;
            }
            SDVertex* newV = new SDVertex(p_new);
            v->child = newV;
            new_vs.push_back(newV);
        }

        // 2. Compute Odd Vertices
        for (auto &shape : sdshapes) {
            for (const SDEdge &e : shape.es) {
                Point p_odd;
                SDVertex *v1 = e.verts[0];
                SDVertex *v2 = e.verts[1];

                if (e.fs[1] == nullptr) { // Boundary edge
                    p_odd = 0.5 * (v1->p + v2->p);
                } else { // Interior edge
                    SDVertex *v3 = e.fs[0]->otherVertex(&e);
                    SDVertex *v4 = e.fs[1]->otherVertex(&e);
                    p_odd = 0.375 * (v1->p + v2->p) + 0.125 * (v3->p + v4->p);
                }
                SDVertex* oddV = new SDVertex(p_odd);
                new_vs.push_back(oddV);

                // Assign oddV to faces
                for(int i=0; i<2; ++i) {
                    SDFace* f = e.fs[i];
                    if(f) {
                        int idx1 = f->getVertIndex(v1);
                        int idx2 = f->getVertIndex(v2);
                        int edgeIdx = -1;
                        if (NEXT(idx1) == idx2) edgeIdx = idx1;
                        else if (NEXT(idx2) == idx1) edgeIdx = idx2;
                        
                        if(edgeIdx != -1) {
                            f->sdverts[edgeIdx] = oddV;
                        }
                    }
                }
            }
        }

        // 3. Re-triangulate
        for (auto &shape : sdshapes) {
            std::vector<SDFace*> new_fs;
            for (SDFace *f : shape.fs) {
                SDVertex *v0 = f->verts[0]->child;
                SDVertex *v1 = f->verts[1]->child;
                SDVertex *v2 = f->verts[2]->child;

                SDVertex *e0 = f->sdverts[0]; // Edge 0-1
                SDVertex *e1 = f->sdverts[1]; // Edge 1-2
                SDVertex *e2 = f->sdverts[2]; // Edge 2-0

                new_fs.push_back(new SDFace(v0, e0, e2));
                new_fs.push_back(new SDFace(e0, v1, e1));
                new_fs.push_back(new SDFace(e1, v2, e2));
                new_fs.push_back(new SDFace(e0, e1, e2));

                delete f; 
            }
            shape.fs = new_fs;
            shape.es.clear(); 
        }

        // 4. Cleanup Old Vertices
        for(auto* v : sdattrib.vs) {
            delete v;
        }
        sdattrib.vs = new_vs;
        
        // Update indices for the new vertices
        for (size_t i = 0; i < sdattrib.vs.size(); ++i) {
            sdattrib.vs[i]->index = i;
        }

        sdattrib.ns.clear(); 
        sdattrib.ts.clear(); 
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