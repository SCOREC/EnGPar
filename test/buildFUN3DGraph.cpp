#include <apfMesh2.h>
#include <cassert>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apf.h>
#include <cstdlib>
#include <stdint.h>
#include <apfNumbering.h>
#include <cstring>
#include <engpar_support.h>

typedef std::vector<apf::MeshEntity*> BL_Verts;

void collapseBoundaryLayer(apf::Mesh*, std::vector<BL_Verts>&);
int numberStacksAndVerts(apf::Mesh*, const std::vector<BL_Verts>&, apf::MeshTag*&);

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  EnGPar_Initialize();
  if ( argc < 6 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <output graph file> <Collapse BL [0/1]> <secondary_dimensions, ...>\n", argv[0]);
    EnGPar_Finalize();
    MPI_Finalize();
    assert(false);
  }

  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);

  bool collapseBL = atoi(argv[4]);
  std::vector<BL_Verts> bl_stacks;
  if (collapseBL) {
    collapseBoundaryLayer(m,bl_stacks);
  }

  printf("Number of stacks found: %lu\n",bl_stacks.size());

  //Create a tag to number vertices
  //All members of a stack have the same number (index of the array)
  apf::MeshTag* vert_ids;
  int num_verts = numberStacksAndVerts(m, bl_stacks, vert_ids);
  
  m->destroyNative();
  apf::destroyMesh(m);

  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");
  EnGPar_Finalize();
  MPI_Finalize();
  return 0;

}

void collapseBoundaryLayer(apf::Mesh* m, std::vector<BL_Verts>& bl_stacks) {
  apf::MeshEntity* v;
  apf::MeshIterator* vitr = m->begin(0);
  //Loop over all vertices
  while ((v = m->iterate(vitr))) {
    apf::ModelEntity* gent = m->toModel(v);
    //Skip any vertex on model elements
    if (m->getModelType(gent)==3)
      continue;
    apf::Adjacent elms;
    m->getAdjacent(v,3,elms);
    int i;
    //Check if the vertex bounds a prism
    for (i = 0; i < elms.size(); ++i) {
      if (m->getType(elms[i]) == apf::Mesh::Type::PRISM)
        break;
    }
    if (i == elms.size())
      continue;

    //Start a BL stack
    BL_Verts stack;
    stack.push_back(v);
    
    apf::Up edges;
    m->getUp(v,edges);
    apf::MeshEntity* e = NULL;
    apf::MeshEntity* next;
    //Proceed up the boundary layer until we can't go any further
    do {
      next = NULL;
      apf::MeshEntity* edge;
      for (int j = 0; j < edges.n; ++j) {
        edge = edges.e[j];
        //Skip the previous edge
        if (edge==e)
          continue;
        apf::Up faces;
        m->getUp(edge,faces);
        int k;
        //Skip any edge that is adjacent to a traingle
        //??? What if a stack is next to a pyramid, but keeps going much further?
        for (k = 0; k < faces.n; ++k) {
          if (m->getType(faces.e[k]) == apf::Mesh::Type::TRIANGLE)
            break;
        }
        if (k != faces.n)
          continue;
        next = edge;
        //Find the new vertex of the stack
        apf::Downward verts;
        m->getDownward(next,0,verts);
        v = (v == verts[0]) ? verts[1] : verts[0];
        stack.push_back(v);
        m->getUp(v,edges);
        break;
      }
      e = next;
    }
    while (e != NULL);

    //Add the stack
    //??? Should we skip stacks of size 1? 2? What do these mean? Are they possible?
    bl_stacks.push_back(stack);
  }
  m->end(vitr);
}


int numberStacksAndVerts(apf::Mesh* m, const std::vector<BL_Verts>& bl_stacks, apf::MeshTag*& vert_ids) {
  vert_ids = m->createIntTag("vertex_ids",1);

  //First mark all bl stack vertices
  for (int i = 0; i < bl_stacks.size(); i++) {
    for (size_t j = 0; j < bl_stacks[i].size(); j++) {
      m->setIntTag(bl_stacks[i][j], vert_ids, &i);
    }
  }

  //Then mark all other vertices
  int id = bl_stacks.size();
  apf::MeshEntity* v;
  apf::MeshIterator* vitr = m -> begin(0);
  while ((v = m -> iterate(vitr))) {
    if (!m->hasTag(v, vert_ids)) {
      m->setIntTag(v, vert_ids, &id);
      id++;
    }
  }
  m->end(vitr);

  printf("There will be %d vertices in the graph\n",id);
  return id;
}
