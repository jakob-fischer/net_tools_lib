/* author: jakob fischer (jakob@automorph.info)
 * date: 9th February 2013
 * description: 
 * Tools for creating different types of complex networks and for
 * analyzing them.
 */

#ifndef NETWORK_TOOLS_H
#define NETWORK_TOOLS_H

#include <vector>
#include <cstdlib>
#include <cmath>
#include <gmp.h>
#include <iml.h>
#include <iostream>
#include <numeric>
#include <algorithm>



class bt_draw_ref {
    size_t accum;
    std::vector<size_t> probs;
	
public:	
	bt_draw_ref() : accum(0) {}
	
	
	bt_draw_ref& add(size_t t) {
	    probs.push_back(t);
	    accum += t;
	    
	    return (*this);	
	}
	
	size_t get_prob(size_t i) {
	    return probs[i];
	}
	
	
	bt_draw_ref& set_prob(size_t id, size_t v) {
	    accum += v;
	    accum -= probs[id];
		
	    probs[id]=v;
	    return (*this);
	}
	
	
	size_t draw() const {
		if(accum == 0) {
		    std::cout << "bt_draw::draw() - accum is zero!" << std::endl;
		    return 0;
		}
			
	    size_t rnd=rand()%accum;
	    size_t z=0;
	    
	    while(probs[z] < rnd) 
	        rnd -= probs[z++];
	
	    return z;
	}
	
	
	size_t size() const {
	    return probs.size();	
	}
	
	
	bt_draw_ref(const std::vector<size_t>& p) : accum(0)  {
	    for(size_t i=0; i<p.size(); ++i)
	        add(p[i]);
	}
};



class bt_draw {
    static const size_t e_max=30; 
	
    size_t last;
    std::vector<size_t> probs[e_max];
    
    
	
public:	
	bt_draw() : last(0) {}
	
	
	
	
	bt_draw& set_prob(size_t id, size_t v) {
	    size_t old=probs[0][id];

	    for(size_t i=0; i<e_max && probs[i].size() != 0; ++i) {
	        probs[i][id/pow(2,i)]  -= old;
	        probs[i][id/pow(2,i)]  += v;	    	    
	    }
	    
	    return (*this);
	}
		
	
	bt_draw& add(size_t t) {
	    probs[0].push_back(0);
	    
	    size_t l=0;
	    while(probs[l].size() != 2*probs[l+1].size() &&
	    	  probs[l].size() != 1) {
	        probs[l+1].push_back(0);
	    	    
	        if(probs[l+1].size() == 1) {
	            ++last;
	            probs[l+1][0] = probs[l][0];
	        }
	        
	        ++l;
	    }
	    
	    // TODO
	    
	    return set_prob(probs[0].size()-1, t);	
	}
	
	size_t get_prob(size_t i) {
	    return probs[0][i];
	}
	

	
	size_t draw() const {
	    if(last == 0) {
	        std::cout << "bt_draw::draw() - last is zero!" << std::endl;
	        return 0;
	    }
			
	    size_t rnd=rand()%probs[last-1][0];
	    
	    size_t z=0;
	    for(size_t i=last-1; i>0; --i) {
	    	// Either not at the right flank or there, but having two childs
	        if(z != probs[i].size()-1 || probs[i].size() == 2*probs[i-1].size()) {
	            if(probs[i-1][z*2] > rnd) {
	                z=z*2;
	            } else {
                        rnd = rnd - probs[i-1][z*2];
	                z=z*2+1;
	            }
	        } else { // Just having one child
	        	z=z*2;
	        }
	    	    
	    }

	    return z;
	}
	
	
	size_t size() const {
	    return probs[0].size();	
	}
	
	
        bt_draw(const std::vector<size_t>& p) : last(0)  {        	
	    for(size_t i=0; i<p.size(); ++i)
	        add(p[i]);
	}
};






/*
 * Helper for "create_pan_sinha".
 * Calculates the size (number of nodes) of the module on level 'level' with the
 * index 'element. It needs a vector 'elem_c' with sizes of elementary modules (level 0)
 * and the parameter 'm' 
 */

size_t bps_module_size(const std::vector<size_t> &elem_c, size_t m,
			      size_t level, size_t element) {
    if(level == 0) 
        return elem_c[element];
    
    if(level == 1) {
        size_t sum=0;
	
	    for(size_t i=0; i<m; ++i) 
	        sum += elem_c[element*m+i];
	
	    return sum;
    }
    
    return bps_module_size(elem_c, m, level-1, 2*element) +
           bps_module_size(elem_c, m, level-1, 2*element+1);  
}


/*
 * Similar to previous function, but calculates the offset. 
 * (between first element of the network and first one in the module)
 */

size_t bps_module_offset(const std::vector<size_t> &elem_c, size_t m,
			      size_t level, size_t element) {
    if(level == 0) {
        size_t sum=0;
	
        for(size_t i=0; i<element; ++i)
            sum += elem_c[i];
	  
        return sum;
    }
	
    if(level == 1)
	    return bps_module_offset(elem_c, m, 0, m*element);
    
    return bps_module_offset(elem_c, m, level-1, 2*element);
}


/*
 * Find the number (0-indexed) of the module on level "level" in which
 * the node "node" is located
 */

size_t bps_module_no(const std::vector<size_t> &elem_c, size_t m,
			      size_t level, size_t node) {
	size_t max_modules=0;
	if(level==0) 
	    max_modules=elem_c.size();
	else 
	    max_modules=elem_c.size()/(m*pow(2,level-1)); 
					  
					  
	for(size_t i=0; i<max_modules; ++i) {
		size_t offset=bps_module_offset(elem_c, m, level, i);
        size_t size=bps_module_size(elem_c, m, level, i);
        
		if(node >= offset  &&
		   node < offset+size)
		    return i; 
	}			  
	
	std::cout << "Error: bps_module_no - did not find module!" << std::endl;
	return 0;
}


/*
 * Helper for 'create_pan_sinha' that checks if two nodes 'n1' and 'n2' are in 
 * the same level-0 module.
 */

bool bps_same_emodule(const std::vector<size_t> & elem_c, size_t n1, size_t n2) {
    size_t first=0;
    size_t sum=0;
    while(sum + elem_c[first] < n1) 
        sum += elem_c[first++];  

    size_t second=0;
    sum=0;
    while(sum + elem_c[second] < n2) 
        sum += elem_c[second++];  
 
    return first==second;  
}


/*
 * Calculates the maximal number of edges for networks with 'nodes' nodes.
 * Depending on if 'self_loop's are allowed and if the network is 'directed'.
 * Multiple edges are not taken into account.
 */

size_t bps_max_edges(size_t nodes, bool self_loop, bool directed) {
    if(self_loop && directed)
        return (nodes*nodes);
    else if(self_loop && !directed)
        return ((nodes*(nodes-1))/2+nodes);
    else if(!self_loop && directed)
        return nodes*nodes-nodes;
    else if(!self_loop && !directed)
        return (nodes*(nodes-1))/2;
}

/*
 * This function returns the (lowest) level on which the two nodes 'a' and 'b'
 * are in the same module. This allows to calculate the probability the 
 * edge betwen them exists.
 */
 
 size_t bps_nodes_level_comp(size_t a, size_t b, 
                             std::vector<size_t>& elem_c, size_t m) {
    size_t lvl=0;
    while(bps_module_no(elem_c, m, lvl, a) !=	
          bps_module_no(elem_c, m, lvl, b))
        ++lvl;							 
					
    return lvl;							 
}


size_t bps_nodes_level_comp_lu(const std::vector<size_t> &modules_no, size_t N,
                               size_t a, size_t b) {
    size_t lvl=0;
    while(modules_no[N*lvl+a] != modules_no[N*lvl+b])
        ++lvl;							 
					
    return lvl;										   
}



/*
 * The function creates a hierarchical-modular network using the algorithm 
 * of pan-sinha. The function allows as well 'multiple' edges, setting of
 * 'self_loop' and 'directed' networks. The model works whith first adding
 * connections inside of modules of the lowest level 1 with probability 'p'. 
 * 'm' of these modules form a module on the level 2. where connections between
 * elements of the lower level modules are done with probability p*r     
 * 
 * This old function distributes the M edges on the different modules. The
 * new one should a global probability distribution on all possible edges and 
 * sample from it.
 */

void create_pan_sinha_old(std::vector< std::pair<size_t, size_t> >& edges, 
                      size_t N, size_t M, size_t h, size_t m, double r, 
                      bool multiple =false, bool self_loop=false, 
                      bool directed=false) {
    size_t n_el_mod=pow(2,h)*m;                   // number of elementary modules 
    size_t n_floor = floor(double(N)/n_el_mod);   // min no nodes in el. modules
    
    std::vector<size_t> elem_c;      // calculate the size of elementary nodes starting by n_floor
    {                                // and adding 1 until N is reached
        size_t N_tmp=N;
   
        for(size_t i=0; i<n_el_mod; ++i) { 
            elem_c.push_back(n_floor);
	    N_tmp -= n_floor;
        }

        for(size_t i=0; i<N_tmp; ++i) 
            ++elem_c[i];
    }
    

    // Calculate maximum number of edges for each level
    std::vector<double> max_edges_level, max_edges_level_r;
    std::vector<size_t> possible_base;    
    
    {   // lowest level
        max_edges_level.push_back(0);
        for(size_t k=0; k<n_el_mod; ++k) {
            size_t nodes=bps_module_size(elem_c, m, 0, k);
            max_edges_level[0] += bps_max_edges(nodes, self_loop, directed);
        }
    }     

    {   // medium level
        max_edges_level.push_back(0);
        for(size_t k=0; k<n_el_mod/m; ++k) {
            size_t nodes=bps_module_size(elem_c, m, 1, k);
            size_t d=bps_max_edges(nodes, self_loop, directed);

            for(size_t j=0; j<m; ++j)
                d -= bps_max_edges(elem_c[m*k+j],self_loop, directed);

            max_edges_level[1] += d;
        }
    }

    // the upper levels
    for(size_t i=1; i<=h; ++i) {
        max_edges_level.push_back(0);
        for(size_t k=0; k<pow(2,h-i); ++k) {
            size_t nodes1=bps_module_size(elem_c, m, i, k*2);
            size_t nodes2=bps_module_size(elem_c, m, i, k*2+1);

            if(directed)
                max_edges_level[i+1] += 2*nodes1*nodes2;
            else 
                max_edges_level[i+1] += nodes1*nodes2;
        }
    }

    // rescale in terms of r
    for(size_t i=0; i<max_edges_level.size(); ++i)
        max_edges_level_r.push_back(max_edges_level[i]*pow(r, i));
    

    // Distribute M edges on all levels according to max_edges_level-distribution
    double sum_mel=std::accumulate(max_edges_level_r.begin(),max_edges_level_r.end(),0);
    size_t M_tmp=M;
    std::vector<size_t> edges_level_partition;
    std::vector<double> edges_level_partition_dif;

    for(size_t i=0; i<max_edges_level.size(); ++i) {
        size_t tmp=floor(M*max_edges_level_r[i]/sum_mel);    
        M_tmp -= tmp;
        edges_level_partition.push_back(tmp);
        edges_level_partition_dif.push_back(M*max_edges_level_r[i]/sum_mel-tmp);
    }

    while(M_tmp > 0) {
        std::vector<double>::iterator it= 
              std::max_element(edges_level_partition_dif.begin(), 
                               edges_level_partition_dif.end());
        size_t next=it-edges_level_partition_dif.begin();
        ++edges_level_partition[next];
        edges_level_partition_dif[next]=0;
        --M_tmp;
    }
       
    std::cout << "Show edge -> level assignment:" << std::endl;
    for(size_t i=0; i<edges_level_partition.size(); ++i)
        std::cout << i << "  :  " << edges_level_partition[i] << "   (" << max_edges_level[i]/pow(r,i) << ")" << std::endl;

    
    // Now generate the edges
    
    {   // lowest level - distribute 
        std::vector<size_t> dist_low_level_edges;
	std::vector<double> dist_low_level_edges_dif;
	size_t elp_tmp=edges_level_partition[0];
	
	for(size_t i=0; i<n_el_mod; ++i) {
	    //std::cout << "elem_c[" << i << "]=" << elem_c[i] << "     bpsm=" << bps_max_edges(elem_c[i], multiple, self_loop, directed) <<  std::endl;
	  
	  
	    double rel=double(bps_max_edges(elem_c[i], self_loop, directed)/max_edges_level[0]*edges_level_partition[0]);
	    dist_low_level_edges.push_back(floor(rel));
	    dist_low_level_edges_dif.push_back(rel-floor(rel));
	    elp_tmp -= floor(rel);
	}
	
	while(elp_tmp > 0) {
            std::vector<double>::iterator it= 
                  std::max_element(dist_low_level_edges_dif.begin(), 
                                   dist_low_level_edges_dif.end());
            size_t next=it-dist_low_level_edges_dif.begin();
            ++dist_low_level_edges[next];
            dist_low_level_edges_dif[next]=0;
            --elp_tmp;
	}
	 
	//  DO IT
	for(size_t i=0; i<n_el_mod; ++i) {
	    size_t offset = bps_module_offset(elem_c, m, 0, i);

	    while(dist_low_level_edges[i] > 0) {
	        size_t n1=offset+rand()%elem_c[i];
	        size_t n2=offset+rand()%elem_c[i];
	        
		if(n1 != n2 || self_loop) {
		    if(directed) {
	                if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n1, n2)) == edges.end() || multiple) {
		                edges.push_back(std::pair<size_t, size_t> (n1, n2) );
		            --dist_low_level_edges[i];
		        }	 
		    } else {
	                if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n1, n2)) == edges.end() &&
			   find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n2, n1)) == edges.end() || multiple) {
		                edges.push_back(std::pair<size_t, size_t> (n1, n2) );
		            --dist_low_level_edges[i];
		        }
		    } 
		}
	    }
	}
    }      
    
    
    {   // medium level - the most complicated ;)
        std::vector<size_t> dist_med_level_edges;
	std::vector<double> dist_med_level_edges_dif;
	size_t elp_tmp=edges_level_partition[1];
	
	for(size_t k=0; k<n_el_mod/m; ++k) {
	    size_t nodes=bps_module_size(elem_c, m, 1, k);
            size_t d=bps_max_edges(nodes, self_loop, directed);

            for(size_t j=0; j<m; ++j)
                d -= bps_max_edges(elem_c[m*k+j] ,self_loop, directed);

	    double rel=double(d/max_edges_level[1]*edges_level_partition[1]);
	    
	    dist_med_level_edges.push_back(floor(rel));
 	    dist_med_level_edges_dif.push_back(rel-floor(rel));
	    elp_tmp -= floor(rel);
	}
	
	while(elp_tmp > 0) {
            std::vector<double>::iterator it= 
                  std::max_element(dist_med_level_edges_dif.begin(), 
                                   dist_med_level_edges_dif.end());
            size_t next=it-dist_med_level_edges_dif.begin();
            ++dist_med_level_edges[next];
            dist_med_level_edges_dif[next]=0;
            --elp_tmp;
	}

	// DO IT
	 
        for(size_t k=0; k<n_el_mod/m; ++k) {
	    size_t offset = bps_module_offset(elem_c, m, 1, k);	
	    size_t nodes=bps_module_size(elem_c, m, 1, k);
	    
	    while(dist_med_level_edges[k] > 0) {
	        size_t n1=offset+rand()%nodes;  
	        size_t n2=offset+rand()%nodes; 
	      
		if(!bps_same_emodule(elem_c, n1, n2)) {
                    if(directed) {
	                if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n1, n2)) == edges.end() || multiple) {
		                edges.push_back(std::pair<size_t, size_t> (n1, n2) );
		            --dist_med_level_edges[k];
		        }	 
		    } else {
	                if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n1, n2)) == edges.end() ||
			   find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n2, n1)) == edges.end() || multiple) {
		                edges.push_back(std::pair<size_t, size_t> (n1, n2) );
		            --dist_med_level_edges[k];
		        }
		    } 		  
		  
		  
		  
		}
	    }
	}
    }
    
    {    // connect the higher levels 
        for(size_t i=1; i<=h; ++i) {
            std::vector<size_t> dist_x_level_edges;
	    std::vector<double> dist_x_level_edges_dif;
	    size_t elp_tmp = edges_level_partition[i+1];
	  
            for(size_t k=0; k<pow(2,h-i); ++k) {
                size_t nodes1=bps_module_size(elem_c, m, i, k*2);
                size_t nodes2=bps_module_size(elem_c, m, i, k*2+1);
                size_t pot_edges=nodes1*nodes2;
                if(directed)
                    pot_edges *= 2;
	
		double rel=double(pot_edges/max_edges_level[i+1]*edges_level_partition[i+1]);
	        dist_x_level_edges.push_back(floor(rel));
	        dist_x_level_edges_dif.push_back(rel-floor(rel));
	        elp_tmp -= floor(rel);
            }
             
            while(elp_tmp > 0) {
                std::vector<double>::iterator it= 
                      std::max_element(dist_x_level_edges_dif.begin(), 
                                       dist_x_level_edges_dif.end());
                size_t next=it-dist_x_level_edges_dif.begin();
                ++dist_x_level_edges[next];
                dist_x_level_edges_dif[next]=0;
                --elp_tmp;
	    }	    
	    
	    // NOW we are doing the connecting
            for(size_t k=0; k<pow(2,h-i); ++k) {
                size_t nodes1=bps_module_size(elem_c, m, i, k*2);
                size_t nodes2=bps_module_size(elem_c, m, i, k*2+1);
		size_t offset=bps_module_offset(elem_c, m, i+1, k);
		
		while(dist_x_level_edges[k] > 0) {
		    size_t n1=offset+rand()%nodes1;
		    size_t n2=offset+nodes1+rand()%nodes2;
		    
		    if(directed) {
		        if(rand()%2) {
			    if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n1, n2)) == edges.end() || multiple) {
		                edges.push_back(std::pair<size_t, size_t> (n1, n2) );
		                --dist_x_level_edges[k];
		            }	  
			} else {
			    if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n2, n1)) == edges.end() || multiple) {
		                edges.push_back(std::pair<size_t, size_t> (n2, n1) );
		                --dist_x_level_edges[k];
		            }	
			}
		    }else {
		        if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n1, n2)) == edges.end() ||
			   multiple) {
		            edges.push_back(std::pair<size_t, size_t> (n1, n2) );
		            --dist_x_level_edges[k];
		        }	  
		    }
		  
		}
            }
        }    
    }
}




/*
 * Creates a 
 * 
 * edges           - reference to edgeliste where links are written to (0-based)
 * N               - number of nodes in network
 * M               - number of edges / links in network
 * allow_multiple  - the same link can occur multiple times
 * self_loop       - self loops are allowed
 * directed        - directed network is created
 * m               - 
 * r               -
 */ 

void create_pan_sinha(std::vector< std::pair<size_t, size_t> >& edges, 
                      size_t N, size_t M, size_t h, size_t m, double r, 
                      bool multiple =false, bool self_loop=false, 
                      bool directed=false) {
    size_t n_el_mod=pow(2,h)*m;                   // number of elementary modules 
    size_t n_floor = floor(double(N)/n_el_mod);   // min no nodes in el. modules
    
    std::vector<size_t> elem_c;      // calculate the size of elementary nodes starting by n_floor

    {                                // and adding 1 until N is reached
        size_t N_tmp=N;
   
        for(size_t i=0; i<n_el_mod; ++i) { 
            elem_c.push_back(n_floor);
	        N_tmp -= n_floor;
        }

        for(size_t i=0; i<N_tmp; ++i) 
            ++elem_c[i];
    }
    
    // Creating lookup table for asociating every node on every level with one module
    std::vector<size_t> module_no;
    size_t module_no_lvl=h+2;
    for(size_t i=0; i<module_no_lvl; ++i) 
		for(size_t j=0; j<N; ++j) 
			module_no.push_back(bps_module_no(elem_c, m, i, j));
    
    
    bt_draw nw_prob;
    for(size_t count=0; count<N*N; ++count) {    
    	    // translate to edge
		size_t i=count/N;
		size_t j=count%N;
		
		if((directed || i <= j) && (self_loop  || i != j))  
			nw_prob.add( pow(r, bps_nodes_level_comp_lu(module_no, N, i, j))*100000 );
	    else 
			nw_prob.add(0);
			
	}
    
    for(size_t a=0; a<M; ++a) {  	
	// draw random number and find asocciated entry
	size_t count=nw_prob.draw();
		
	// translate to edge
	size_t i=count/N;
	size_t j=count%N;
			
        edges.push_back(std::pair<size_t, size_t> (i, j) );
            
        if(!multiple)  // if no multiples allowed set prob. for next run to 0
        	nw_prob.set_prob(count, 0);	
	}
}


/*
 * Creates a simple barabasi-albert-network with N nodes and M edges. Edges are written
 * to the vector given as first parameter. Edges are indexed starting with zero.
 * 
 * edges           - reference to edgeliste where links are written to (0-based)
 * N               - number of nodes in network
 * M               - number of edges / links in network
 * allow_multiple  - the same link can occur multiple times
 * self_loop       - self loops are allowed
 * directed        - directed network is created
 * 
 * 
 */ 
 
void create_barabasi_albert(std::vector< std::pair<size_t, size_t> > &edges, 
                            size_t N, size_t M, bool multiple=false, 
                            bool self_loop=false, bool directed=true) {
    // vector dynamically containing the functionality (degree centrality) 
    // of all nodes
    std::vector<size_t> fun;
    size_t fun_sum=0;             // sum over fun (for normalization)
    for(size_t i=0; i<N; ++i)
        fun.push_back(0);
  
    // Implement first edge between node 0 and 1
    edges.push_back( std::pair<size_t, size_t>(0,1) );
    ++fun[0];
    ++fun[1];
    fun_sum += 2;
    M -= 1;
    
    // Proceeding with node 2 and higher
    for(size_t i=2; i<N; ++i) {
		// - distribute number of edges that remain to be placed on number
		// of nodes that have to be "created"
		// - maximum number of edges connected to current node is number of
		// already existing nodes
        size_t no_n_edges = M/(N-i) > i ? i : M/(N-i);
	
	    for(size_t k=0; k<no_n_edges; ++k) {
	        // try creating edge number k (until done is true)
	        bool done=false;
	            
	        do {
				// select node randomly with probability of their functionality
	            size_t l=0;
	            long r=rand()%fun_sum;
		        
		        while(r >= fun[l]) {
		            r -= fun[l];
		            ++l;
		        }
	      
	            // proceed if no self loop or self loop allowed
	            if(l != i || self_loop) {
	                if(directed) {
	                    if(find(edges.begin(), edges.end(), 
			                    std::pair<size_t, size_t>(i, l)) == edges.end() ||
		                   multiple) {
							   
							// directed edges are directed randomly   
			                if(rand()%2)
		                        edges.push_back(std::pair<size_t, size_t> (i, l) );
			                else
			                    edges.push_back(std::pair<size_t, size_t> (l, i) );
			    
			                ++fun[i];
			                ++fun[l];
			                fun_sum += 2;
			    
		                    done = true;
		                }		  
		            } else {
	                    if(find(edges.begin(), edges.end(), 
			               std::pair<size_t, size_t>(i, l)) == edges.end() ||
		                   multiple) {
	                        // for non directed the higher numbered node is always
	                        // put in first
		                    edges.push_back(std::pair<size_t, size_t> (i, l) );
			                ++fun[i];
			                ++fun[l];
			                fun_sum += 2;
		                    done = true;
		                }
		            }
                }
	        } while(!done);
	    }
	
	    M -= no_n_edges;
    }
}




/*
 * Creates a simple watts-strogatz-network with N nodes and M edges. Edges are written
 * to the vector given as first parameter. Edges are indexed starting with zero.
 * 
 * edges           - reference to edgeliste where links are written to (0-based)
 * N               - number of nodes in network
 * M               - number of edges / links in network
 * allow_multiple  - the same link can occur multiple times
 * self_loop       - self loops are allowed
 * directed        - directed network is created
 */

void create_watts_strogatz(std::vector< std::pair<size_t, size_t> > &edges, 
                           size_t N, size_t M, double beta, bool multiple=false, 
                           bool self_loop=false, bool directed=true) {
    // the parameters self_loop, directed and multiple are only checked
    // for the shuffling
    size_t M_tmp=M;     // counts the remaining number of links to be created
    
    for(size_t dif=1; M_tmp > 0; ++dif) {
        for(size_t base=0; base<N && M_tmp > 0; ++base) {
            size_t n1=base;
	        size_t n2=(base+dif)%N;
	  
	        // Random direction if directed, if not there is no harm done
	        if(rand()%2) {
	             edges.push_back( std::pair<size_t, size_t> (n1, n2) );  
	        } else {
	            edges.push_back( std::pair<size_t, size_t> (n2, n1) );
	        }
	    
	        --M_tmp;
	    }   
    }
 
    // Shuffle second node of edge 'i' with probability beta
    for(size_t i=0; i<edges.size(); ++i) {
        if(double(rand())/RAND_MAX < beta) {
	        size_t first=edges[i].first;
	        bool found=false;
	    
	        // draw random numbers for second node until you find one 
	        // that matches the constraints
	        do {
	            size_t second=rand()%N;
		
		        if(second != first || self_loop) {
		            if(directed) {
		                if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(first, second)) == edges.end() ||
		                   multiple) {
		                    edges[i].second = second;
		                    found = true;
			            }		      
		            } else {
		                if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(first, second)) == edges.end() &&
		                   find(edges.begin(), edges.end(), std::pair<size_t, size_t>(second, first)) == edges.end() ||
		                   multiple) {
		                    edges[i].second = second;
		                    found = true;
			            }
	                }
		        }
	        } while (!found);
	    }
    }
}



/*
 * Creates a simple erdos-renyi-network with N nodes and M edges. Edges are written
 * to the vector given as first parameter. Edges are indexed starting with zero.
 * 
 * edges           - reference to edgeliste where links are written to (0-based)
 * N               - number of nodes in network
 * M               - number of edges / links in network
 * allow_multiple  - the same link can occur multiple times
 * self_loop       - self loops are allowed
 * directed        - directed network is created
 */

void create_erdos_renyi(std::vector< std::pair<size_t, size_t> > &edges, 
                        size_t N, size_t M, bool allow_multiple=false, 
                        bool self_loop=false, bool directed = false) {
    for(size_t j=edges.size(); j<M; ++j) {
        // Found a link that is allowed given the existing edges?
        bool found=false;
	
	    do {  
	        size_t first=rand()%N;  
            size_t second=rand()%N;
		
	        // If applicable check for self loops and multiple occorance
	        if(first != second || self_loop) {
	            if(allow_multiple) { // If multiple links is allowed no reason to check for direction
	                found = true;	 
		            edges.push_back(std::pair<size_t, size_t> (first, second) );
	            } else if(directed) { // If directed and no multiples allowed check once
                    if(find(edges.begin(), edges.end(), 
		                std::pair<size_t, size_t>(first, second)) == edges.end()) {
		                edges.push_back(std::pair<size_t, size_t> (first, second) );
		                found = true;
		            }	
                } else {              // Else link can have two forms
                    if(find(edges.begin(), edges.end(), 
		               std::pair<size_t, size_t>(first, second)) == edges.end() &&
		               find(edges.begin(), edges.end(), 
		               std::pair<size_t, size_t>(second, first)) == edges.end()) {
		                edges.push_back(std::pair<size_t, size_t> (first, second) );
		                found = true;
		            }		  
		        }  
	        }
	    } while (!found);
    }
}



/*
 * For an Erdos-Renyi-Network this function this function selects edges to 
 * couple edges. This can be used when a mixed 1 -> 1 / 2->2 - reaction
 * network should be created with an Erdos-Renyi statistics of the reactants 
 * graph. There are just 'C' pairs of links selected (without replacement).
 *
 * couples    - reference to a vector of pairs to that the selected couples are
                written to 0-based 
 * C          - number of couples that are drawn
 * edges      - pairs vector containing the original networks edges
 * 
 * TODO: Should this algorithm be extended to allow maintaining the network
 *       constraints (self_loop, ...) for the sustrate graph
 */

void couple_erdos_renyi(std::vector< std::pair<size_t, size_t> > &couples,
                        size_t C,
                        std::vector< std::pair<size_t, size_t> > &edges) {
    std::vector<bool> taken_f;  // which links are already part of a link pair
    for(size_t i=0; i<edges.size(); ++i)
        taken_f.push_back(false);

    while(C > 0) {  // C pairs have to be found
        size_t first=rand()%edges.size();	
	    size_t second=rand()%edges.size();	
	
        // both partners have to be different and free...    
	    if(first != second && !taken_f[first] && !taken_f[second]) {
	        couples.push_back(std::pair<size_t, size_t> (first, second));
	        taken_f[first]=true;
	        taken_f[second]=true;
	        --C;	
	    }
    }
}


/*
 * For an Watts-Strogatz-Network this function this function selects edges to 
 * couple edges. This can be used when a mixed 1 -> 1 / 2->2 - reaction
 * network should be created with an Watts-Strogatz statistics of the reactants 
 * graph. There are just 'C' pairs of links selected (without replacement).
 *
 * couples    - reference to a vector of pairs to that the selected couples are
                written to 0-based 
 * C          - number of couples that are drawn
 * edges      - pairs vector containing the original networks edges
 * N          - Parameter used for generating the ws-model (needed to maintain clustering)
 * beta       - shuffling parameter in watz-strogatz-model necessary so newly added
 *              links (in imaginary substrate graph) have beta as shuffling probability          
 * 
 * TODO: measure beta from network? Check if algorithm always terminates! 
 * TODO: Should this algorithm be extended to allow maintaining the network
 *       constraints (self_loop, ...) for the sustrate graph    
 */
 
void couple_watts_strogatz(std::vector< std::pair<size_t, size_t> > &couples,
                           size_t C, std::vector< std::pair<size_t, size_t> > &edges,
                           size_t N, double beta) {
    size_t M=edges.size();
    size_t close_v = ceil(double(M)/N); // needed to count which additional "connections" are far

    std::vector<bool> taken_f;  // which links are already part of a link pair
    for(size_t i=0; i<M; ++i)
        taken_f.push_back(false);

    size_t far_links = C*beta;  // TODO: Maybe measure beta by hand?

    while(C > 0) {
        size_t i=rand()%M;
        size_t j=rand()%M;
        if(i != j && !taken_f[i] && !taken_f[j]) {
            size_t count=0;         // number of far edges introduced by this couple 
            if(edges[i].first-edges[j].second > close_v)
                ++count;

            if(edges[j].first-edges[i].second > close_v)
                ++count;

            if(count <= far_links) {   // Dont have to many far links
                if(count == 0 && C < far_links) {   // Avoid adding to few far links
	                couples.push_back(std::pair<size_t, size_t> (i, j));
	                taken_f[i]=true;
	                taken_f[j]=true;

				    C -= 1;
				    far_links -= count;
				}
			}
  
        }
    }
}


/*
 * For a Barabasi-Albert-Network this function this function selects edges to 
 * couple edges. This can be used when a mixed 1 -> 1 / 2->2 - reaction
 * network should be created with an Watts-Strogatz statistics of the reactants 
 * graph. There are just 'C' pairs of links selected (without replacement).
 *
 * couples    - reference to a vector of pairs to that the selected couples are
                written to 0-based 
 * C          - number of couples that are drawn
 * edges      - pairs vector containing the original networks edges
 * N          - Parameter used for generating the ws-model (needed to maintain clustering)
 * beta       - shuffling parameter in watz-strogatz-model necessary so newly added
 *              links (in imaginary substrate graph) have beta as shuffling probability          
 * 
 * TODO: Should this algorithm be extended to allow maintaining the network
 *       constraints (self_loop, ...) for the sustrate graph
 *       also: should the change of functionality through virtual links
 *             be considered in probability of forming future couples?
 * TODO: How should the probability of adding a certain couple be estimated.
 *       maximum product of "new edges" functionality or sum of products? 
 */
 
void couple_barabasi_albert(std::vector< std::pair<size_t, size_t> > &couples,
                            size_t C, std::vector< std::pair<size_t, size_t> > &edges,
                            size_t N) {
    std::vector<size_t> fun_n;
    std::vector<bool> edge_available;
    size_t fun_sum=0;
    size_t M=edges.size();
    for(size_t i=0; i<N; ++i)
        fun_n.push_back(0);
        
    for(size_t j=0; j<M; ++j) {
	    ++fun_n[edges[j].first];	
		++fun_n[edges[j].second];
		fun_sum += 2;
		edge_available.push_back(true);
	}
		
    std::vector<size_t> nw_prob_ac;
    for(size_t k=0; k<M*M; ++k)
	    nw_prob_ac.push_back(0);	
    
    for(size_t run=0; run<C; ++run) {
		// setting up cummulative probability distribution
		size_t acc=0;
        for(size_t c=0; c<nw_prob_ac.size(); ++c) {
			size_t i=c/M;
			size_t j=c%M;
			    
			if(i != j && edge_available[i] && edge_available[j]) {
			    size_t e1_f = edges[i].first;	
			    size_t e1_s = edges[i].second;	
			    size_t e2_f = edges[j].first;	
			    size_t e2_s = edges[j].second;
			    acc += (fun_n[e1_f]+fun_n[e1_s]+fun_n[e2_f]+fun_n[e2_s]);
			} 
			                
			nw_prob_ac[c]=acc;
	    }	
	    
	    // draw random number and find asocciated entry
		size_t sum=nw_prob_ac[nw_prob_ac.size()-1];
		double rnd=rand()%sum;
			
		size_t count=0;
		while(nw_prob_ac[count] < rnd)
			++count;
			    
		// translate to edge
		size_t i=count/N;
		size_t j=count%N;
			
        couples.push_back(std::pair<size_t, size_t> (i, j) ); 
        edge_available[i]=false;
        edge_available[j]=false;
        //TODO change functionality of nodes here?           
	}							 
}


/*
 * 
 * 
 * 
 */

void couple_pan_sinha(std::vector< std::pair<size_t, size_t> > &couples,
                      size_t C, std::vector< std::pair<size_t, size_t> >& edges, 
                      size_t N, size_t h, size_t m, double r) {
	size_t n_el_mod=pow(2,h)*m;                   // number of elementary modules 
    size_t n_floor = floor(double(N)/n_el_mod);   // min no nodes in el. modules
    
    std::vector<size_t> elem_c;      // calculate the size of elementary nodes starting by n_floor

    {                                // and adding 1 until N is reached
        size_t N_tmp=N;
   
        for(size_t i=0; i<n_el_mod; ++i) { 
            elem_c.push_back(n_floor);
	        N_tmp -= n_floor;
        }

        for(size_t i=0; i<N_tmp; ++i) 
            ++elem_c[i];
    }					  


   // Creating lookup table for asociating every node on every level with one module
    std::vector<size_t> module_no;
    size_t module_no_lvl=h+2;
    for(size_t i=0; i<module_no_lvl; ++i) 
		for(size_t j=0; j<N; ++j) 
			module_no.push_back(bps_module_no(elem_c, m, i, j));
						  
						  
    std::vector<bool> edge_available;
    size_t M=edges.size();
        
    for(size_t j=0; j<M; ++j) 
		edge_available.push_back(true);
		
    bt_draw nw_prob;
    for(size_t k=0; k<M*M; ++k)
	    nw_prob.add(0);	
    
    for(size_t run=0; run<C; ++run) {
		// setting up cummulative probability distribution
		size_t acc=0;
		
        for(size_t c=0; c<nw_prob.size(); ++c) {
			size_t i=c/M;
			size_t j=c%M;
			    
			if(i != j && edge_available[i] && edge_available[j]) {
			    size_t e1_f = edges[i].first;	
			    size_t e1_s = edges[i].second;	
			    size_t e2_f = edges[j].first;	
			    size_t e2_s = edges[j].second;
			    nw_prob.set_prob(c, ((pow(r, bps_nodes_level_comp_lu(module_no, N, e1_f, e2_s)) +
			             pow(r, bps_nodes_level_comp_lu(module_no, N, e2_f, e1_s)))*100000));
			} 
	    }	
	    
	    // draw random number and find asocciated entry
		size_t count=nw_prob.draw();
		
		// translate to edge
		size_t i=count/N;
		size_t j=count%N;
			
        couples.push_back(std::pair<size_t, size_t> (i, j) ); 
        edge_available[i]=false;
        edge_available[j]=false;        
	}							    								 						 						  
}





/*
 * Creates an Erdos Renyi network with an stochastical algorithm with N nodes
 * and on average M edges. From M and the maximal number of edges a global 
 * probability for creating edges is calculated.
 */

void create_erdos_renyi_stoch(std::vector< std::pair<size_t, size_t> > &edges, 
                              size_t N, size_t M, bool self_loop=false, 
                              bool directed = false) {
    // For each case, calculate first the probability for a arbitrary link
    // of existing d, so that the expectation value for number of links is M
    // and then iterate over all possible links

    if(!self_loop && !directed) {
        double d=double(M)/(N*(N-1)*0.5);  
      
        for(size_t i=0; i<N; ++i) 
	        for(size_t j=i+1; j<N; ++j)
	            if(double(rand())/RAND_MAX < d)
		            edges.push_back(std::pair<size_t, size_t> (i, j) );        
    } else if (self_loop && !directed) {
        double d=double(M)/(N*(N-1)*0.5+N);  
      
        for(size_t i=0; i<N; ++i) 
	        for(size_t j=i; j<N; ++j)
	            if(double(rand())/RAND_MAX < d)
		            edges.push_back(std::pair<size_t, size_t> (i, j) );        
    } else if (!self_loop && directed) {
        double d=double(M)/(N*N-N); 
      
        for(size_t i=0; i<N; ++i) 
	        for(size_t j=0; j<N; ++j)
	            if(double(rand())/RAND_MAX < d && i!=j)
		            edges.push_back(std::pair<size_t, size_t> (i, j) );      
    } else if (self_loop && directed) {
        double d=double(M)/(N*N);  
      
        for(size_t i=0; i<N; ++i) 
	        for(size_t j=0; j<N; ++j)
	            if(double(rand())/RAND_MAX < d)
		            edges.push_back(std::pair<size_t, size_t> (i, j) );
    }
}





/* bad written class for simple bidirectional network analysis. not 
 * developed and propperly commented at the moment because R and igraph 
 * can do better.
 */


class bd_network {
    size_t nodes, edges;
    bool* ad_m;
    size_t* sh_path;
    
    std::vector<size_t> degrees, deg_hist, sh_path_hist, cluster_id, cluster_size;
    std::vector<double> deg_dist, sh_path_dist, clustering_coeff;
    size_t min_deg, max_deg, max_sh_path, smallest_cluster, biggest_cluster;
    double avg_deg, avg_sh_path, biggest_cluster_ratio, avg_clustering_coeff;
    
  
public:
    /*
     * INPUT OUTPUT CODE
     *
     */
    void add_edge(size_t a, size_t b) {
        if(ad_m[a*nodes+b])
            return;
	  
	++edges;
        ad_m[a*nodes+b]=true;
	ad_m[b*nodes+a]=true;
    }
    
  
    bd_network(size_t nd) : nodes(nd), edges(0), avg_deg(0), avg_sh_path(0), 
	    biggest_cluster_ratio(0), min_deg(0), max_deg(0), max_sh_path(0), 
	    smallest_cluster(0), biggest_cluster(0), avg_clustering_coeff(0) {
        ad_m = new bool[nodes*nodes];
	sh_path = new size_t[nodes*nodes];
	
	for(size_t i=0; i<nodes*nodes; ++i) {
	    ad_m[i]=false;
	    sh_path[i]=0;    
	}
    }
    
    
    bd_network(size_t N, std::vector< std::pair<size_t, size_t> > &edg) : nodes(N), 
            edges(0), avg_deg(0), avg_sh_path(0), biggest_cluster_ratio(0), min_deg(0), 
            max_deg(0), max_sh_path(0), smallest_cluster(0), biggest_cluster(0), avg_clustering_coeff(0)  {
        ad_m = new bool[nodes*nodes];
	sh_path = new size_t[nodes*nodes];
	
	for(size_t i=0; i<nodes*nodes; ++i) {
	    ad_m[i]=false;
	    sh_path[i]=0;    
	}           
      
        for(size_t i=0; i<edg.size(); ++i) 
	    add_edge(edg[i].first, edg[i].second);
    }
    
    
    ~bd_network() {
        delete[] ad_m; 
	delete[] sh_path;
    } 
    
    
    bool is_connected() {
        for(size_t i=0; i<nodes*nodes; ++i)
	    if(sh_path[i] > edges)
	        return false;
	    
	return true;
    }
  
    std::vector<size_t>& get_degrees() {  return degrees;  }
    
    std::vector<size_t>& get_degree_hist() {  return deg_hist;  }
    
    std::vector<double>& get_degree_dist() {  return deg_dist;  }
    
    size_t get_shortest_path(size_t a, size_t b) {  return sh_path[b*nodes+a];  }
    
    std::vector<size_t>& get_path_hist() {   return sh_path_hist;  }
    
    std::vector<double>& get_path_dist()  {  return sh_path_dist;  }
    
    size_t get_min_degree() {  return min_deg;  }
    
    size_t get_max_degree() {  return max_deg;  }
    
    double get_avg_degree() {  return avg_deg;  }
    
    size_t get_max_path() {  return max_sh_path;  }
    
    size_t get_diameter() {  return max_sh_path;  }
    
    double get_avg_path() {  return avg_sh_path;  }
    
    double get_biggest_cluster_ratio() {  return biggest_cluster_ratio;  }
    
    size_t get_smallest_cluster() {  return smallest_cluster;  }
    
    size_t get_biggest_cluster() {  return biggest_cluster;  }
    
    size_t get_cluster_number() {  return cluster_id.size(); }
        
    std::vector<double>& get_clustering_coeff() {  return clustering_coeff;  } 
    
    double get_avg_clustering_coeff() {  return avg_clustering_coeff;  }
    
  
    /*
     * ANALYSIS CODE
     *
     */
    
    
    double get_steady_state_A(size_t a, size_t b, std::vector<double>& concentration) {
        long i, j, bd, s;
        mpz_t *mp_B, *mp_N;

        /* generate a n x m random left hand side matrix A */
        size_t n = nodes;
        size_t m = nodes+1;
        long* A = new long[n*m];    

	// initializing
	for(size_t i=0; i<n*m; ++i)
	    A[i] = 0;
	
	// Creating stoichiometry matrix for (inner) network interaction
	for(size_t i=0; i<nodes; ++i) {
	    for(size_t j=i+1; j<nodes; ++j) {
	        if(ad_m[i*nodes+j]) {	  
	            A[i*m+j] += 1;
	            A[i*m+i] -= 1;
	            A[j*m+i] += 1;
	            A[j*m+j] -= 1;
		}
	    }
	}
	
	// Adding boundary condition / flow
	A[a*m+m-1]=1;
        A[b*m+m-1]=-1;
	
	

        s = nullspaceLong (n, m, A, &mp_N);
        

    
	mpz_t C_11, C_12, C_21, C_22, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, det;
	mpz_init(C_11);   mpz_init(C_12);   mpz_init(C_21);   mpz_init(C_22);   mpz_init(tmp1); 
	mpz_init(tmp2);   mpz_init(tmp3);   mpz_init(tmp4);   mpz_init(tmp5);   mpz_init(tmp6);   mpz_init(det);
	
	mpz_set(C_11, mp_N[a * s + 0]);
	mpz_set(C_12, mp_N[a * s + 1]);
	mpz_set(C_21, mp_N[b * s + 0]);
	mpz_set(C_22, mp_N[b * s + 1]);
	
        delete [] A;
  

	mpz_mul(tmp1, C_11, C_22);
	mpz_mul(tmp2, C_21, C_12);
	mpz_sub(det, tmp1, tmp2);
	
	
	mpz_set(tmp3, mp_N[nodes * s + 0]);
	mpz_set(tmp4, mp_N[nodes * s + 1]);
	mpz_mul(tmp1, C_22, tmp3);
	mpz_mul(tmp2, C_21, tmp4);
	mpz_sub(tmp3, tmp2, tmp1);
	
	mpq_t v;
	mpq_init(v);
        mpq_set_num(v, tmp3);
	mpq_set_den(v, det);
	double rs=mpq_get_d(v);
	
	{   // Calculate the solution vector
	    concentration.clear();
	    
	    for(size_t i=0; i<=nodes; ++i) {
	        double d=0.0;
	      	
		mpz_set(tmp3, mp_N[i * s + 0]);
	        mpz_set(tmp4, mp_N[i * s + 1]);
	        mpz_mul(tmp1, C_22, tmp3);
	        mpz_mul(tmp2, C_21, tmp4);
	        mpz_sub(tmp3, tmp2, tmp1);
	
                mpq_set_num(v, tmp3);
	        mpq_set_den(v, det);
	        d=mpq_get_d(v);
		
	        concentration.push_back(d);
	    }
	}
	
        for (i = 0; i < m * s; i++)
            mpz_clear (mp_N[i]);
        free (mp_N);
	
	// mpz cleanup
	mpz_clear(C_11);   mpz_clear(C_12);   mpz_clear(C_21);   mpz_clear(C_22);
	mpz_clear(tmp1);   mpz_clear(tmp2);   mpz_clear(tmp3);   mpz_clear(tmp4);	
	mpz_clear(tmp5);   mpz_clear(tmp6);   mpz_clear(det);    mpq_clear(v);
	
	return rs;
    }
    
    
   
    
    double get_steady_state_B(size_t a, size_t b, std::vector<double>& concentration) {
        double alpha=get_steady_state_A(a, b, concentration);
        
	for(size_t i=0; i<concentration.size(); ++i)
	    concentration[i] *= (alpha)/(1+2*alpha);
	  
      
	return (alpha)/(1+2*alpha);
    }
    
  
    void calculate_degrees() {  
        degrees.clear();
	deg_hist.clear();
	deg_dist.clear();
      
        // Calculate degree of every node
        for(size_t i=0; i<nodes; ++i) {
	    degrees.push_back(0);  
	    
	    for(size_t j=0; j<nodes; ++j)
	        if(ad_m[i*nodes+j])
		    ++degrees[i];
	}      
	
	// Calculate minimal, maximal and average degree
	avg_deg=0;
	min_deg=max_deg=degrees[0];
	
	for(size_t i=0; i<nodes; ++i) {
	    avg_deg += degrees[i];  
	    
	    if(degrees[i] < min_deg) min_deg=degrees[i];
	    if(degrees[i] > max_deg) max_deg=degrees[i];
	}
	
	avg_deg /= nodes;
	
	// Calculate degree distribution (first histogramm than norm. dist.)	
	for(size_t i=0; i<=max_deg; ++i) {
	    deg_hist.push_back(0);
	    
	    for(size_t j=0; j<nodes; ++j) 
	        if(degrees[j] == i)
		    ++deg_hist[i];
		
	    deg_dist.push_back(double(deg_hist[i])/nodes);
	}
      
    }
  
 
    /*
     * Calculating the shortest paths between all pairs of network nodes...
     */

    void calculate_sh_paths() {
        sh_path_dist.clear();
	sh_path_hist.clear();
      
        // Calculate length of paths	
        for(size_t i=0; i<nodes; ++i) {
	    // The row is inicialiced
            for(size_t j=0; j<nodes; ++j) 
	        sh_path[i*nodes+j]=nodes+2;  // Marks a not yet known distance
	    
	    sh_path[i*nodes+i]=0;  // initial condition
	    
	    // now the distance is calculated recursively
	    //calculate_dmap(i, i);
	    
	    size_t current=0;
	    size_t highest=nodes+2;
	    while(highest==nodes+2) {
	        highest = 0;
		
		//std::cout << "current is " << current << std::endl;
		
		for(size_t j=0; j<nodes; ++j) {
		    if(sh_path[i*nodes+j] == current) {		      
		        for(size_t k=0; k<nodes; ++k) 
			    if(ad_m[j*nodes+k] && sh_path[i*nodes+k] > current)
			        sh_path[i*nodes+k] = current+1;
		    } else if(sh_path[i*nodes+j] > highest) {
		        highest = sh_path[i*nodes+j];
		    }
		}
	        
	        ++current;
	    }
	}
		
	
	// Calculate shortest path length distribution / average / maximum 
	max_sh_path=0;
	size_t count=0;
	
	// Iterating all pairs of nodes
        for(size_t i=0; i<nodes; ++i) {
            for(size_t j=0; j<nodes; ++j) {
	        if(sh_path[i*nodes+j] != 0 &&        // Path length has to be finite
		   sh_path[i*nodes+j] <= edges) {
		   
		   avg_sh_path += sh_path[i*nodes+j];
		   ++count;
		
		   // updating maximum shortest path
		   if(max_sh_path < sh_path[i*nodes+j])
		       max_sh_path = sh_path[i*nodes+j];
		}
	    }
	} 
	
	// Initializing path histogramm vector
	for(size_t i=0; i<=max_sh_path; ++i)
	    sh_path_hist.push_back(0);
	
	// Calculating histogramm
	for(size_t i=0; i<nodes; ++i) 
            for(size_t j=0; j<nodes; ++j) 
	        if(sh_path[i*nodes+j] != 0 &&
		   sh_path[i*nodes+j] <= edges) 
		    ++sh_path_hist[sh_path[i*nodes+j]];  
	
	// .. average path length..
	avg_sh_path /= count;
	
	// ..and shortest path length distribution.
	for(size_t i=0; i<=max_sh_path; ++i)
	    sh_path_dist.push_back(double(sh_path_hist[i])/count);
    }
  
  
    /*
     * Method calculates the clustering coefficient for each node and the average
     * clustering coefficient.
     *
     * Degrees have to be calculated first!
     */
     
    void calculate_clustering_coeff() {
        clustering_coeff.clear();
        avg_clustering_coeff=0;
	
	
	for(size_t i=0; i<nodes; ++i)
	    clustering_coeff.push_back(0);
	    
	for(size_t i=0; i<nodes; ++i) 
	    for(size_t j=0; j<nodes; ++j) 
	        if(ad_m[i*nodes+j]) 
		    // Found pair of connected nodes. Now look for nodes that are 
		    // connected to both....  
		    for(size_t k=0; k<nodes; ++k)
		        if(ad_m[k*nodes+i] && ad_m[k*nodes+j]) 
			    clustering_coeff[k] += 1.0;
	
	// Normalizing the counted size...
	for(size_t i=0; i<nodes; ++i) 
	    if(degrees[i] > 1)
	        clustering_coeff[i] /= (degrees[i]*(degrees[i]-1));
	    else 
	        clustering_coeff[i] = -1;
	
	
	// calculate average clustering coefficient
	avg_clustering_coeff=0;
	size_t cc_count=0;
	for(size_t i=0; i<clustering_coeff.size(); ++i)
	    if(clustering_coeff[i] >= 0) { 
	        avg_clustering_coeff += clustering_coeff[i];
	        ++cc_count;
	    }
	
	avg_clustering_coeff /= cc_count;
    }
  
  
  
    void m_draw_clusters_rec(size_t node, size_t next_color=0) {
        // next_color==0 bedeutet es sollen alle cluster bestimmt werden...
        if(next_color==0) {
            for(size_t i=0; i<nodes; ++i)
	        if(cluster_id[i] == 0)
	            m_draw_clusters_rec(i, ++next_color);
      
	     return;
        }
  
        // Ist der aktuelle Node schon gefärbt wird die Rekursion abgebrochen
        if(cluster_id[node] != 0)
            return;
  
        // Setze die Farbe des aktuellen Nodes (node) auf next_color
        cluster_id[node]=next_color;
  
        // färbe rekursiv alle mit node verbundenen nodes
        for(size_t i=0; i<nodes; ++i) 
	  if(ad_m[node*nodes+i])
	        m_draw_clusters_rec(i, next_color);
      
    }
  
    
    void calculate_con_clusters() {
        cluster_id.clear();
	cluster_size.clear();
	
	for(size_t i=0; i<nodes; ++i)
	    cluster_id.push_back(0);
      
	m_draw_clusters_rec(0, 0);
	
	// calculate number of clusters
	size_t no_clusters=0;
	for(size_t i=0; i<nodes; ++i)
	    if(cluster_id[i] > no_clusters)
	        no_clusters = cluster_id[i];
	
	// calculating cluster sizes
	for(size_t i=0; i<no_clusters; ++i)
	    cluster_size.push_back(0);
	
	for(size_t i=0; i<nodes; ++i) 
	    ++cluster_size[cluster_id[i]-1];
	
	// calculating maximum and minimum cluster size
        smallest_cluster=edges;
	biggest_cluster=0;
	
	for(size_t i=0; i<cluster_size.size(); ++i) {
	    if(cluster_size[i]<smallest_cluster)
	        smallest_cluster=cluster_size[i];
	    
	    if(cluster_size[i]>biggest_cluster)
	        biggest_cluster=cluster_size[i]; 
	}
	
	biggest_cluster_ratio = double(biggest_cluster)/nodes; 
    }
  
  
    void do_analyis() {
        std::cout << "deg..";  	std::cout.flush();
	calculate_degrees();
        std::cout << "sp..";	std::cout.flush();
	calculate_sh_paths();
        std::cout << "cc..";	std::cout.flush();
	calculate_clustering_coeff();
        std::cout << "con..";	std::cout.flush();
	calculate_con_clusters();
	std::cout << std::endl;
     } 
};





















#endif // NETWORK_TOOLS_H
