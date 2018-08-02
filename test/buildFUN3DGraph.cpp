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
#include <ngraph.h>

typedef std::vector<apf::MeshEntity*> BL_Verts;

void collapseBoundaryLayer(apf::Mesh*, std::vector<BL_Verts>&);
int numberStacksAndVerts(apf::Mesh*, const std::vector<BL_Verts>&, apf::MeshTag*&);
void gatherGraphVertices(apf::Mesh*, const std::vector<BL_Verts>&, apf::MeshTag*,
                         agi::lid_t, agi::gid_t*&, agi::wgt_t*&);
agi::lid_t gatherGraphEdges(apf::Mesh*, const std::vector<BL_Verts>&, apf::MeshTag*, int dim,
                            agi::gid_t*& edge_ids, agi::lid_t*& degs, agi::gid_t*& pins,
                            agi::wgt_t*&);

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

  //Construct the vertices of the graph
  agi::gid_t* verts;
  agi::wgt_t* vert_weights;
  gatherGraphVertices(m, bl_stacks, vert_ids, num_verts, verts, vert_weights);
  agi::Ngraph* g = agi::createEmptyGraph();
  g->constructVerts(true, num_verts, verts, vert_weights);
  delete [] verts;
  delete [] vert_weights;

  //Construct the edges of the graph
  agi::gid_t* edge_ids;
  agi::lid_t* degs;
  agi::gid_t* pins;
  agi::wgt_t* edge_weights;
  int dimension = 1;
  agi::lid_t num_edges = gatherGraphEdges(m, bl_stacks, vert_ids, dimension, edge_ids, degs,
                                          pins, edge_weights);
  g->constructEdges(num_edges, edge_ids, degs, pins, edge_weights);
  delete [] edge_ids;
  delete [] degs;
  delete [] pins;
  delete [] edge_weights;
  
  //TODO: construct ghosts (for now only worrying about serial)

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

void gatherGraphVertices(apf::Mesh* m, const std::vector<BL_Verts>& bl_stacks, apf::MeshTag* vert_ids, agi::lid_t nv, agi::gid_t*& verts, agi::wgt_t*& wgts) {

  verts = new agi::gid_t[nv];
  wgts = new agi::wgt_t[nv];

  size_t i;
  for (i =0; i < bl_stacks.size(); i++) {
    verts[i] = i;
    wgts[i] = bl_stacks[i].size();
  }
  for (; i < nv; i++) {
    verts[i] = i;
    wgts[i] = 1;
  }
}


typedef std::set<int> HyperEdge;
//A hash functor for the a set of ints. Simply add the hashes of the entries.
//Note: A better hash function may improve results if this approach has too many collisions
class HEHash {
public:
  size_t operator()(const HyperEdge& e) const {
    int h=0;
    std::hash<int> hasher;
    for(HyperEdge::iterator itr = e.begin(); itr != e.end(); itr++) {
      h+=hasher(*itr);
    }
    return h;
  }
};

//A simple operator equals for the sets for use in the unordered map
inline bool operator==(const HyperEdge& e1, const HyperEdge& e2) {
  if (e1.size()!=e2.size())
    return false;
  //HyperEdge are ordered sets so if they are equal the elements will be in the same order
  HyperEdge::iterator itr1;
  HyperEdge::iterator itr2;
  for (itr1 = e1.begin(), itr2 = e2.begin(); itr1 != e1.end(); ++itr1, ++itr2) {
    if (*itr1 != *itr2)
      return false;
  }
  return true;
}


agi::lid_t gatherGraphEdges(apf::Mesh* m, const std::vector<BL_Verts>& bl_stacks,
                            apf::MeshTag* vert_ids, int dim, agi::gid_t*& edge_ids,
                            agi::lid_t*& degs, agi::gid_t*& pins, agi::wgt_t*& wgts) {
  //An unordered map to store hyperedges (sets of vertex ids) temporarily to count duplicates
  typedef std::unordered_map<HyperEdge, int, HEHash> HEMap;
  HEMap hyperedges;

  agi::gid_t num_pins = 0;
    
  apf::MeshEntity* bridge;
  apf::MeshIterator* itr = m->begin(dim);
  while ((bridge = m->iterate(itr))) {
    apf::Downward verts;
    int nd = m->getDownward(bridge,0,verts);
    HyperEdge e;
    for (int i = 0; i < nd; ++i) {
      int id;
      m->getIntTag(verts[i],vert_ids,&id);
      e.insert(id);
    }
    std::pair<HEMap::iterator, bool> insert_pair  = hyperedges.insert(std::make_pair(e,1));
    if (!insert_pair.second)
      insert_pair.first->second++;
    else {
      num_pins+=e.size();
    }
  }
  m->end(itr);

  agi::lid_t num_edges = hyperedges.size();
  printf("There will be %d edges and %ld pins in the graph\n", num_edges, num_pins);

  edge_ids = new agi::gid_t[num_edges];
  degs = new agi::lid_t[num_edges];
  pins = new agi::gid_t[num_pins];
  wgts = new agi::wgt_t[num_edges];

  HEMap::iterator heItr;
  int index = 0;
  int pin_index = 0;
  for (heItr = hyperedges.begin(); heItr != hyperedges.end(); ++heItr) {
    edge_ids[index] = index;
    degs[index] = heItr->first.size();
    HyperEdge::iterator hItr;
    for (hItr = heItr->first.begin(); hItr != heItr->first.end(); ++hItr) {
      pins[pin_index++] = *hItr;
    }
    wgts[index++] = heItr->second;
  }

  return num_edges;
}
