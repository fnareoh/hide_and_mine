//
//  shorten_Z.cpp
//  shorten
//
//  Created by 陈慧娉 on 30/01/2019.
//  Copyright © 2019 陈慧娉. All rights reserved.
//

#include "shorten_Z.hpp"
#include <string>
#include <cstring>
#include <cstdlib>
#include <time.h>
#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <stack>  
#include <map>
#include <typeinfo>
//#include "/usr/local/Cellar/igraph/0.7.1_3/include/igraph/igraph.h"
#include "/usr/local/include/igraph/igraph.h"

using namespace std;

void graph_print(igraph_t graph){
    igraph_eit_t eit;
    igraph_integer_t from, to;
    igraph_real_t w;
    igraph_integer_t eid;
    
    printf("\n BEGIN***graph_print**(node ids start from 0)************ \n");
    igraph_eit_create(&graph,  igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
    
    while (!IGRAPH_EIT_END(eit))
    {
        eid = IGRAPH_EIT_GET(eit);
        
        
        igraph_edge(&graph, eid,&from, &to);
        
        printf("[%d](%d, %d) \n", eid, from, to);
        
        
        IGRAPH_EIT_NEXT(eit);
    }
    printf("END*********************** \n");
    
    igraph_eit_destroy(&eit);
    
}
void print_A(unordered_map<int, unordered_set<int> >& A)
{
    cout<<"Array A: "<<endl;
    for(unordered_map<int, unordered_set<int> >::const_iterator it=A.begin();it!=A.end();++it)
    {
        cout<<"Letter "<<it->first<<"-> paths: ";
        for(unordered_set<int>::const_iterator it2=(it->second).begin();it2!=(it->second).end();++it2)
        {
            cout<<*it2<<" ";
        }
        cout<<endl;
    }
    
}
void update_A(unordered_map<int, unordered_set<int> >& A, vector<int> path, int path_id)
{
    for(vector<int>::iterator it=path.begin();it!=path.end();++it)
    {
        unordered_map<int, unordered_set<int> >::iterator it2=A.find(*it);
        
        
        if ( it2 == A.end() )
        {
            unordered_set<int> new_list;
            new_list.insert(path_id);
            A[*it]=new_list;
        }
        else
        {
            it2->second.insert(path_id);
        }
    }
}

void print_set_of_paths(vector<vector<int> >& set_of_paths_P)
{
    int cnt=0;
    for(vector<vector<int> >::iterator it=set_of_paths_P.begin();it!=set_of_paths_P.end();++it)
    {
        cout<<"p["<<cnt<<"]: ";
        for(vector<int>::iterator it2=(*it).begin();it2!=(*it).end();++it2)
        {
            cout<<*it2<<" ";
        }
        cout<<endl;
        ++cnt;
    }
}
void join_path_and_circle(unordered_map<int, unordered_set<int> >& A, vector<vector<int> >& set_of_paths_P, int node, vector<int>& intersecting_path, 
                          int intersecting_path_id, vector<int>& circle)
{
    vector<int> joined_path;
    
    vector<int> tmp;
    
    vector<int>::iterator node_position_in_path, node_position_in_circle;
    
    for(vector<int>::iterator it=intersecting_path.begin();it!=intersecting_path.end();++it)
    {
        if(*it!=node)
        {
            joined_path.push_back(*it);
        }
        else
        {
            
            joined_path.push_back(*it);
            node_position_in_path=it;
            break;
        }
        
    }
    
    for(vector<int>::iterator it=circle.begin();it!=circle.end();++it)
    {
        if(*it==node)
        {
            node_position_in_circle=it;
            break;
        }
    }
    
    ++node_position_in_circle;
    
    if(node_position_in_circle!=circle.end())
    {
        for(vector<int>::iterator it=node_position_in_circle;it!=circle.end();++it)
        {
            joined_path.push_back(*it);
        }
        
    }
    
    --node_position_in_circle;
    
    for(vector<int>::iterator it=circle.begin();it!=node_position_in_circle;++it)
    {
        joined_path.push_back(*it);
    }
    joined_path.push_back(*node_position_in_circle);
    
    ++node_position_in_path;
    if(node_position_in_path!=intersecting_path.end())
    {
        for(vector<int>::iterator it=node_position_in_path;it!=intersecting_path.end();++it)
        {
            joined_path.push_back(*it);
        }
    }
    
    set_of_paths_P.at(intersecting_path_id)=joined_path;
    
    update_A(A, joined_path, intersecting_path_id);
}
void find_a_path_with_node(unordered_map<int, unordered_set<int> >& A, vector<vector<int> >& set_of_paths_P, int node, vector<int>& path, int& intersecting_path_id)
{
    unordered_map<int, unordered_set<int> >::iterator got = A.find(node);
    
    if(got!=A.end() && !(got->second).empty())
    {
        intersecting_path_id=*(got->second.begin());
        path=set_of_paths_P.at(intersecting_path_id);
    }
}

void create_A(unordered_map<int, unordered_set<int> >& A, vector<vector<int> >& set_of_paths_P)
{
    int cnt=0;
    for(vector<vector<int> >::iterator it=set_of_paths_P.begin();it!=set_of_paths_P.end();++it)
    {
        for(vector<int>::iterator it2=(*it).begin();it2!=(*it).end();++it2)
        {
            
            unordered_map<int, unordered_set<int> >::iterator got = A.find(*it2);
            
            if ( got == A.end() )
            {
                unordered_set<int> new_list;
                new_list.insert(cnt);
                A[*it2]=new_list;
            }
            else
            {
                got->second.insert(cnt);
            }
            
        }
        ++cnt;
    }
}

void vector_print(igraph_vector_t *v) {
    long int i;
    for (i=0; i<igraph_vector_size(v); i++) {
        printf(" %li", (long int) VECTOR(*v)[i]);
    }
    printf("\n");
}

void delete_cycle(igraph_t& g, vector<int>& cycle, int last_node)
{
    vector<int>::const_iterator it;
    
    for(it=cycle.begin();it!=cycle.end();++it)
    {
        int start_node=*it;
        ++it;
        if(it==cycle.end())
        {
            break;
        }
        
        int end_node=*it;
        
        igraph_es_t es;
        igraph_vector_t v;
        igraph_vector_init(&v,2);
        igraph_integer_t from=start_node, to=end_node;
        igraph_vector_set(&v,0,from);
        igraph_vector_set(&v,1,to);
        
        igraph_es_pairs(&es, &v, true);
        
        igraph_delete_edges(&g, es);
        
        igraph_vector_destroy(&v);
        
        --it;
    }
    
    --it;
    
    igraph_es_t es;
    igraph_vector_t v;
    igraph_vector_init(&v,2);
    igraph_integer_t from=*it, to=last_node;
    igraph_vector_set(&v,0,from);
    igraph_vector_set(&v,1,to);
    
    igraph_es_pairs(&es, &v, true);
    
    igraph_delete_edges(&g, es);
    
    igraph_vector_destroy(&v);
    
}

unsigned int shorten(int l, char *y, vector<string> &Z)
{
    string my_y(y);
    
    char **B=new char*[my_y.length()+1];
    
    char *temp;
    int i=0;
    
    temp=strtok(y,"#");
    
    if(temp!=NULL)
    {
        string temp_str(temp);
        
        
        B[0]=new char[temp_str.length()];
        B[0] =temp;
    }
    
    while (temp!=NULL){
        i++;
        temp = strtok(NULL,"#");
        if(temp!=NULL)
        {
            string temp_str2(temp);
            B[i]=new char[temp_str2.length()];
            B[i] = temp;
        }
    }
    
    
    
    string t;
    
    set<string> B_ps;
    
    
    multimap<string, int> pref_string;
    
    
    multimap<string, int> suff_string;
    
    for (int m=0;m<i;m++)
    {
        t=B[m];
        
        B_ps.insert(t.substr(0,l));
        B_ps.insert(t.substr(t.length()-l,l));
        
        pref_string.insert(pair<string,int>(t.substr(0,l),m));
        suff_string.insert(pair<string,int>(t.substr(t.length()-l,l),m));
    }
    
    
    vector<pair<int, int> > B_prime;
    string s;
    
    for (int m=0;m<i;m++)
    {
        set<string>::iterator it;
        t = B[m];
        int counter_pre=  1;
        int counter_suf = 1;
        int resualt1 = 0;
        int resualt2 = 0;
        for(it = B_ps.begin ();it!= B_ps.end ();it++){
            s = *it;
            if(!resualt1){
                if(t.substr(0,l).compare(s)==0){
                    resualt1 = counter_pre;
                }
                else{
                    counter_pre++;
                }
            }
            if(!resualt2){
                
                if(t.substr(t.length()-l,l).compare(s)==0){
                    
                    resualt2 = counter_suf;
                }
                else{
                    counter_suf++;
                }
            }
        }
        
        B_prime.push_back(make_pair(resualt1, resualt2));
    }
    
    int no_graph_nodes=0;
    
    for (int i=0;i<B_prime.size();i++)
    {
        if(B_prime[i].first>no_graph_nodes || B_prime[i].second>no_graph_nodes)
            no_graph_nodes=max(B_prime[i].first,B_prime[i].second);
    }
    
    
    igraph_t g;
    igraph_vector_t v1;
    int ret;
    
    igraph_i_set_attribute_table(&igraph_cattribute_table);
    
    
    igraph_vector_init(&v1, B_prime.size()*2);
    
    int cnt=0;
    for(vector<pair<int,int> >::const_iterator it=B_prime.begin();it!=B_prime.end();++it)
    {
        VECTOR(v1)[cnt]=(int)(*it).first-1;
        cnt++;
        
        VECTOR(v1)[cnt]=(int)(*it).second-1;
        cnt++;
    }
    
    igraph_create(&g, &v1, 0, 1);
    
    
    igraph_es_t es;
    
    igraph_es_all(&es,IGRAPH_EDGEORDER_ID);
    
    igraph_eit_t eit;
    
    igraph_eit_create(&g, es, &eit);
    
    
    vector<int> in_less_out;
    
    igraph_vs_t vs;
    igraph_vs_all(&vs);
    
    igraph_vit_t vit;
    igraph_integer_t size;
    
    int err=igraph_vit_create(&g, vs, &vit);
    if(err!=0)
        cout<<"vit error: "<<err<<endl;
    
    while(!IGRAPH_VIT_END(vit)) {
        
        int vertex_id=( int)IGRAPH_VIT_GET(vit);
        igraph_vector_t result_indegree, result_outdegree;
        igraph_vector_init(&result_indegree, 1);
        igraph_vector_init(&result_outdegree, 1);
        
        igraph_vs_t vid;
        igraph_vs_1(&vid, vertex_id);
        
        igraph_degree(&g, &result_indegree, vid,IGRAPH_IN,1);
        igraph_degree(&g, &result_outdegree, vid, IGRAPH_OUT,1);
        
        if(VECTOR(result_indegree)[0]<VECTOR(result_outdegree)[0])
            in_less_out.push_back(vertex_id);
        
        igraph_vector_destroy(&result_indegree);
        igraph_vector_destroy(&result_outdegree);
        IGRAPH_VIT_NEXT(vit);
    }
    
    vector<vector<int> > set_of_paths_P;
    
    srand (time(NULL));
    
    bool started_path=false;
    int random_node=-1;
    
    while(!in_less_out.empty())
    {
        vector<int> path_p;
        
        random_node=in_less_out.at(rand()%in_less_out.size());
        
        path_p.push_back(random_node);
        
        igraph_vector_t out_neigh_of_random;
        igraph_vector_init(&out_neigh_of_random, 1);
        
        
        igraph_neighbors(&g, &out_neigh_of_random, random_node, IGRAPH_OUT);
        
        int out_degree=igraph_vector_size(&out_neigh_of_random);
        
        
        if(out_degree>0)
        {
            int random_out_neighbor_index=rand()%igraph_vector_size(&out_neigh_of_random);
            int random_out_neighbor=VECTOR(out_neigh_of_random)[random_out_neighbor_index];
            
            
            igraph_es_t es;
            igraph_vector_t v;
            igraph_vector_init(&v,2);
            igraph_integer_t from=random_node, to=random_out_neighbor;
            igraph_vector_set(&v,0,from);
            igraph_vector_set(&v,1,to);
            
            igraph_es_pairs(&es, &v, true);
            
            igraph_delete_edges(&g, es);
            
            igraph_vector_destroy(&v);
            
            path_p.push_back(random_out_neighbor);
            
            
            igraph_vector_t new_outneigh;
            igraph_vector_init(&new_outneigh, 1);
            
            while(true)
            {
                igraph_neighbors(&g, &new_outneigh, random_out_neighbor, IGRAPH_OUT);
                int out_degree=igraph_vector_size(&new_outneigh);
                
                if(out_degree==0)
                {
                    
                    set_of_paths_P.push_back(path_p);
                    break;
                }
                else
                {
                    
                    int random_out_neighbor_index2=rand()%igraph_vector_size(&new_outneigh);
                    int random_out_neighbor2=VECTOR(new_outneigh)[random_out_neighbor_index2];
                    
                    
                    igraph_es_t es;
                    igraph_vector_t v;
                    igraph_vector_init(&v,2);
                    igraph_integer_t from=random_out_neighbor, to=random_out_neighbor2;
                    
                    igraph_vector_set(&v,0,from);
                    igraph_vector_set(&v,1,to);
                    
                    igraph_es_pairs(&es, &v, true);
                    
                    igraph_delete_edges(&g, es);
                    igraph_vector_destroy(&v);
                    
                    path_p.push_back(random_out_neighbor2);
                    
                    
                    random_out_neighbor=random_out_neighbor2;
                }
            }
        }
        else
        {
            set_of_paths_P.push_back(path_p);
        }
        
        igraph_vector_destroy(&out_neigh_of_random);
        
        
        igraph_vector_t result_indegree, result_outdegree;
        igraph_vector_init(&result_indegree, 1);
        igraph_vector_init(&result_outdegree, 1);
        
        igraph_vs_t vid;
        igraph_vs_1(&vid, random_node);
        
        igraph_degree(&g, &result_indegree, vid,IGRAPH_IN,1);
        igraph_degree(&g, &result_outdegree, vid, IGRAPH_OUT,1);
        
        if(VECTOR(result_indegree)[0]>=VECTOR(result_outdegree)[0])
        {
            std::vector<int>::iterator position = std::find(in_less_out.begin(), in_less_out.end(), random_node);
            if (position != in_less_out.end())
            {
                in_less_out.erase(position);
            }
            
        }
    }
    
    unordered_map<int, unordered_set<int> > A;
    
    create_A(A, set_of_paths_P);
    
    while(igraph_ecount(&g)!=0)
    {
        
        int random_node=rand()%igraph_vcount(&g);
        
        
        igraph_vs_t vid;
        igraph_vs_1(&vid, random_node);
        
        igraph_vector_t out_neigh_of_random;
        igraph_vector_init(&out_neigh_of_random, 1);
        
        
        int current_id=random_node;
        
        unordered_set<int> visited;
        visited.insert(current_id);
        
        vector<int> current_path;
        current_path.push_back(current_id);
        bool cycle_found=false;
        
        do
        {
            igraph_neighbors(&g, &out_neigh_of_random, current_id, IGRAPH_OUT);
            int out_degree=igraph_vector_size(&out_neigh_of_random);
            
            if(out_degree==0)
            {
                break;
            }
            
            int random_out_neighbor_index=rand()%igraph_vector_size(&out_neigh_of_random);
            int random_out_neighbor_id=VECTOR(out_neigh_of_random)[random_out_neighbor_index];
            
            unordered_set<int>::const_iterator visited_it=visited.find(random_out_neighbor_id);
            
            if(visited_it!=visited.end())
            {
                
                vector<int>::const_iterator first;
                
                for(vector<int>::const_iterator it=current_path.begin();it!=current_path.end();++it)
                {
                    if(*it==random_out_neighbor_id)
                    {
                        first=it;
                    }
                }
                
                current_path.erase (current_path.begin(),first);
                
                delete_cycle(g,current_path,random_out_neighbor_id);
                
                
                cycle_found=true;
                break;
                
            }
            else
            {
                visited.insert(random_out_neighbor_id);
                current_path.push_back(random_out_neighbor_id);
            }
            
            current_id=random_out_neighbor_id;
            
        }while(current_id!=random_node);
        
        vector<int> intersecting_path;
        int intersecting_path_id=-1;
        
        if(cycle_found)
        {
            for(vector<int>::const_iterator it=current_path.begin();it!=current_path.end();++it)
            {
                int current_circle_node=*it;
                
                find_a_path_with_node(A, set_of_paths_P, current_circle_node, intersecting_path, intersecting_path_id);
                
                if(!intersecting_path.empty())
                {
                    
                    join_path_and_circle(A, set_of_paths_P, current_circle_node, intersecting_path, intersecting_path_id, current_path);
                    
                    break;
                }
                
            }
            
            if(intersecting_path_id==-1)
            {
                
                current_path.push_back(current_path.front());
                
                set_of_paths_P.push_back(current_path);
                
                int current_path_id=set_of_paths_P.size()-1;
                
                update_A(A, current_path, current_path_id);
                
            }
        }
        
    }
    
    
    int path_id=0;
    int last_path_id=set_of_paths_P.size()-1;
    
    set<int> used_indices_of_B_prime;
    
    for(vector<vector<int> >::const_iterator it=set_of_paths_P.begin();it!=set_of_paths_P.end();++it)
    {
        
        vector<string> B_ps_vec(B_ps.begin(),B_ps.end());
        
        int last=(*it).size()-1;
        int cnt=0;
        
        
        bool deleted_string=false;
        
        int one_before_last=(*it).size()-1;
        
        string string_for_current_prefix;
        
        for(vector<int>::const_iterator it2=(*it).begin();it2!=(*it).end();++it2)
        {
            if(cnt<one_before_last)
            {
                vector<int>::const_iterator it_new;
                if(it2!=(*it).end())
                    it_new=it2+1;
                
                int B_prime_index=0;
                for(vector<pair<int,int> >::const_iterator it3=B_prime.begin();it3!=B_prime.end();++it3)
                {
                    
                    set<int>::iterator used_it=used_indices_of_B_prime.find(B_prime_index);
                    
                    if(it3->first== ((*it2)+1) && (it3->second)==((*it_new)+1) && used_it==used_indices_of_B_prime.end())
                    {
                        
                        string current=B[B_prime_index];
                        if(cnt==0)
                        {
                            Z.push_back(current);
                            if(cnt==one_before_last-1)
                                Z.push_back("#");
                        }
                        else
                        {
                            Z.push_back(current.substr(l,current.length()-l));
                            if(cnt==one_before_last-1)
                                Z.push_back("#");
                        }
                        
                        used_indices_of_B_prime.insert(B_prime_index);
                        break;
                    }
                    B_prime_index++;
                }
            }
            cnt++;
        }
        
        
        path_id++;
    }
    return 1;
}


