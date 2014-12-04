#include<iostream>
#include<cstdlib>
#include<ctime>
#include<cmath>

#define C1 1.5
#define C2 1.1 
using namespace std;


typedef struct tag3{       //to be used as a node for the tree
        struct tag3* parent;
        int size;        
}treenode;

struct arrelem{   //an element for point array
        double x;
        double y;
        int X;
        int Y;    
};

typedef struct tag{         //a node of adjacency list
        struct tag *next;
        int index;
        double d;        
}node;

typedef struct tag2{        //an edge element
        struct tag2 *next;
        int u;
        int v;
        double cost;        
}edge;

typedef struct tag4{         //an element ;to be used for hash table 
	struct tag4 *next;
	int index;
}indexnode;

typedef node *adjl;

double dist(arrelem arr[],int a,int b) //computes distances between points
{
       return sqrt((arr[a].x-arr[b].x)*(arr[a].x-arr[b].x)+(arr[a].y-arr[b].y)*(arr[a].y-arr[b].y));       
}

void heapify(edge* edgearr,int i,int n)  //heapifies array
{
     int l,r,min;
     edge t;
     
     while(1)
     {
             l=2*i+1;r=2*i+2;
             
             if(l>=n) return;
             
             min=((r==n)||(edgearr[l].cost<edgearr[r].cost))?l:r;
             
             if(edgearr[min].cost>edgearr[i].cost) return;
             
             t=edgearr[i];
             edgearr[i]=edgearr[min];
             edgearr[min]=t;
             i=min;        
     }
}

void makeheap(edge* edgearr,int size)     //converts array to heap
{
     int i;
     
     for(i=size/2;i>=0;i--)
          heapify(edgearr,i,size);
}

void deletemin(edge* edgearr,int size)  //deletes minimum element of heap
{
     edgearr[0]=edgearr[size-1];
     heapify(edgearr,0,size-1);
}

typedef treenode* Tree;

Tree kruskal(edge* edgearr,int n,int size)     //implementation of kruskal algorithm
{
    Tree header=(Tree)malloc(n*sizeof(treenode));  //header initialized
    edge e;
    int i,flag;
    Tree t1,t2,s;
    
    for(i=0;i<n;i++)
    {
         header[i].size=0                        ;
         header[i].parent=NULL;
    }    

    while(size>=1)  //while more edges remain
    {
         e=edgearr[0];  //minimum cost edge
         t1=&header[e.u];
         t2=&header[e.v];
         
         while(t1->parent!=NULL)  //go to root
             t1=t1->parent;
         
         while(t2->parent!=NULL)  //go to root
             t2=t2->parent;
         
         if(t1!=t2)
         {
              if(t1->size<t2->size)
              {     
		      t1=&header[e.u];
                    while(t1->parent!=NULL)  //path compression by making the root of the v the parent of all parents of u and u itself 
		            {		    
		               s=t1->parent;
		               t1->parent=t2;
	                   t1=s;
	                }	
              (t2->size)+=((t1->size)+1);      //size of v incremented
                    t1->parent=t2;
                    t1->size=0;                     
              }
              else
              {
		      t2=&header[e.v];
                    while(t2->parent!=NULL)      //path compression by making the root of the v the parent of all parents of u and u itself 
		            {		    
                        s=t2->parent;
		                t2->parent=t1;
	                    t2=s;
                    }
              (t1->size)+=((t2->size)+1);     //size of u incremented
                    t2->parent=t1;
                    t2->size=0;
              }
        }
    deletemin(edgearr,size);  //to get next minimum
    size--;
    }    
    return header;     
}

indexnode** inithash(int m)
{
    int i,j;
    indexnode **hash=(indexnode**)malloc((m)*sizeof(indexnode*));//hash memory allocated
	
    for(i=0;i<m;i++)
	hash[i]=(indexnode*)malloc((m)*sizeof(indexnode));	
	
    for(i=0;i<m;i++)                       //hash initialized
                    for(j=0;j<m;j++)
	                {
                        hash[i][j].next=NULL;
	                    hash[i][j].index=-1;
                    }
    return hash;            
}

void initializepoints(arrelem* arr,indexnode** hash,int m,int n)
{   
    double x,y;
    indexnode *tempnewnode,*newnode;
    int i;
    int X,Y;
    arrelem p;
    
    for(i=0;i<n;i++)  //randomly generated points allocated
    {
	                   x=(double)rand()/RAND_MAX;
	                   y=(double)rand()/RAND_MAX;
	                   X=(int)(x*m);
	                   Y=(int)(y*m);
	                   
                       p.x=x;
	                   p.y=y;
	                   p.X=X;
	                   p.Y=Y;
	                   arr[i]=p;  //point array formation
              		   
                       newnode=(indexnode*)malloc(sizeof(indexnode));   
			           newnode->index=i;
	                   tempnewnode=hash[X][Y].next;
			           hash[X][Y].next=newnode;   //hash formation
			           newnode->next=tempnewnode;
    }         
}

adjl initlist(int n)  //adjacency list initialization
{
    int i;
    adjl ad=(adjl)malloc(n*sizeof(node));
    
    for(i=0;i<n;i++)
    {
         (ad[i]).next=NULL;
         (ad[i]).index=-1;
    }
    return ad;     
}

int formlist(adjl ad,arrelem* arr,indexnode** hash,int n,int m,double d1)  //formation of the contents of the adjacency list
{
    indexnode *l;
    adjl nnode,temp;
    int  cnt=0;
    int i,j,k;
    
    for(i=0;i<n;i++)
    {
                for(j=arr[i].X-(int)ceil(m*d1);j<=arr[i].X+(int)ceil(m*d1);j++)//limits
                         for(k=arr[i].Y-(int)ceil(m*d1);k<=arr[i].Y+(int)ceil(m*d1);k++)
                                if(j>=0&&j<m&&k>=0&&k<m)//to keep in bounds
                                        for(l=hash[j][k].next;l!=NULL;l=l->next)  //hash traversal
                                                if(dist(arr,i,l->index)<=d1&&(i!=l->index)&&((arr[l->index].x-arr[i].x>=0&&arr[l->index].y-arr[i].y>=0)||(arr[l->index].x-arr[i].x>=0&&arr[l->index].y-arr[i].y<=0)))
                                                {
                                                     nnode=(adjl)malloc(sizeof(node));
                                                     nnode->index=l->index;
                                                     nnode->d=dist(arr,i,l->index);
                                                     temp=ad[i].next;
                                                     ad[i].next=nnode;
                                                     nnode->next=temp;   
                                                     cnt++;                                      
                                                }                              
                
    }
    return cnt;
}

struct stagepara{
       int sum;
       int num;       
};

stagepara stageparameters(Tree stage1,int n)  //returns various stage parameters
{
    stagepara res;
    res.sum=0;res.num=0;
    int i;
    
    for(i=0;i<n;i++)
         if(stage1[i].size!=0)
         {
            res.sum+=stage1[i].size;
            res.num++;                     
         }
    return res;       
}

struct formnewverticessettype{
       int* newvertices;
       int newsize;       
};

int formlist2(adjl ada,arrelem* arr,int* newvertices,indexnode** hash,int newsize,int m,double d1,double d2)  //forms list for stage 2
{
    int cnt,i,j,k;
    indexnode* l;
    adjl nnode,temp;
    cnt=0;
    for(i=0;i<newsize;i++)
    {
                for(j=arr[newvertices[i]].X-(int)ceil(m*d2);j<=arr[newvertices[i]].X+(int)ceil(m*d2);j++) //limits
                      for(k=arr[newvertices[i]].Y-(int)ceil(m*d2);k<=arr[newvertices[i]].Y+(int)ceil(m*d2);k++)
                             if(j>=0&&j<m&&k>=0&&k<m)
                                    for(l=hash[j][k].next;l!=NULL;l=l->next)//hash traversal
                                           if(dist(arr,newvertices[i],l->index)<=d2&&dist(arr,newvertices[i],l->index)>d1&&(newvertices[i]!=l->index)&&((arr[l->index].x-arr[newvertices[i]].x>=0&&arr[l->index].y-arr[newvertices[i]].y>=0)||(arr[l->index].x-arr[newvertices[i]].x>=0&&arr[l->index].y-arr[newvertices[i]].y<=0)))  //restrictions
                                           {
                                                 nnode=(adjl)malloc(sizeof(node));
                                                 nnode->index=l->index;
                                                 nnode->d=dist(arr,newvertices[i],l->index);
                                                 temp=ada[newvertices[i]].next;
                                                 ada[newvertices[i]].next=nnode;
                                                 nnode->next=temp;   
                                                 cnt++;                                          
                                           }                              
    }    
    return cnt;
}

formnewverticessettype formnewverticesset(Tree stage1,int n)      //forms the set of new vertices for next stage
{
    formnewverticessettype newset;
    int max=0,maxi,newsize=0;
    int i;
    Tree t1,t2;
    
    for(i=0;i<n;i++)  //to find the giant tree
    {
        if(stage1[i].size>max)
        {
        max=stage1[i].size;
        maxi=i;
        }                
    }
    Tree tmax=&stage1[maxi];
    for(i=0;i<n;i++)  //to get newsize
    {
        t1=&stage1[i];
        while(t1->parent!=NULL)
        t1=t1->parent;
        
        if(t1!=tmax)  
        newsize++;     
    }

    int* newvertices=(int*)malloc(newsize*sizeof(int));  //newvertices memory allocated
    int j=0;
    for(i=0;i<n;i++)
    {
        t1=&stage1[i];
        
        while(t1->parent!=NULL)
        t1=t1->parent;
        
        if(t1!=tmax)  
        {
              newvertices[j]=i;j++;
        } 
    }                 
    newset.newvertices=newvertices;
    newset.newsize=newsize;
    return newset;
}

edge* formnewedgearr(adjl ada,int* newvertices,edge* edgearrcopy,int size,int cnt,int newsize)  //forms ne set of edges
{
    int i,j;
    adjl temp;
    edge* newedgearr=(edge*)malloc((cnt+size)*sizeof(edge));
    for(i=0;i<size;i++)
    newedgearr[i]=edgearrcopy[i];

    free(edgearrcopy);
    
    j=size;
    for(i=0;i<newsize;i++)  //new edge array formed
    {
                   
        for(temp=ada[newvertices[i]].next;temp!=NULL;temp=temp->next)
        {
        newedgearr[j].u=newvertices[i];
        newedgearr[j].v=temp->index;
        newedgearr[j].cost=temp->d;                                
        j++;
        }
    }      
    return newedgearr;
}


void rajasekaran(arrelem* arr,indexnode** hash,long int n,int m)
{
    int i,j,cnt;
    double d1=C1/sqrt(n),d;      //first radius obtained

    adjl ad=initlist(n);    //list initialized
    
    cnt=formlist(ad,arr,hash,n,m,d1);   //list formed

    adjl ptr;

    edge* edgearr=(edge*)malloc(cnt*sizeof(edge));  //edge array initialized
	
    j=0;
    for(i=0;i<n;i++)   //edgearr  formed
    {   
       for(ptr=ad[i].next;ptr!=NULL;ptr=ptr->next)
       {
           edgearr[j].u=i;
           edgearr[j].v=ptr->index;
           edgearr[j].cost=ptr->d;
           j++;                                        
       }             
    }
    
    makeheap(edgearr,cnt);  //edge array heapified
    
    edge* edgearrcopy=(edge*)malloc(cnt*sizeof(edge));  //to be used for new edge array
    
    for(i=0;i<cnt;i++)   //edge array copy made
    {
        edgearrcopy[i]=edgearr[i];
    }
    
    Tree stage1=kruskal(edgearr,n,cnt);    //stage 1 kruskal applied
    
    stagepara stageres=stageparameters(stage1,n);       //stage 1 parameters

    cout<<"\n+++ First phase of EMST...\nNo of edges in G is "<<cnt<<"\nNo of edges in T is "<<stageres.sum<<"\nNumber of trees = "<<stageres.num;
    
    double d2=C2*log(n)/sqrt(n);  //second radius found
    
    int size=cnt;
    adjl ada=initlist(n);   //list initialized
    
    formnewverticessettype newset;
    
    newset=formnewverticesset(stage1,n);  //new set of vertices parameters found
    int* newvertices=newset.newvertices;
    int newsize=newset.newsize;
    
    cnt=formlist2(ada,arr,newvertices,hash,newsize,m,d1,d2);  //number of edges from adjacency list determined; list formed
    
    edge* newedgearr=formnewedgearr(ada,newvertices,edgearrcopy,size,cnt,newsize);  //newedge array formed
    
    int totsize=cnt+size;
    makeheap(newedgearr,totsize);      //heap of new edge array made

    Tree stage2=kruskal(newedgearr,n,totsize);
    int sum=0,num=0;int max=0;
    
    for(i=0;i<n;i++)          //largest tree size determined
    {
        if(stage1[i].size>max)
        {
        max=stage1[i].size;
        }                
    }
    
   stageres=stageparameters(stage2,n); //final output obtained for stage 2
	
    cout<<"\n\n+++ Second phase of EMST...\nLargest tree is of size:"<<max<<"\nNumber of edges in G is:"<<totsize<<"\nNumber of edges in T:"<<stageres.sum<<"\nNumber of trees:"<<stageres.num<<"\n";
         
}

int main()  //main function
{
    long int n;
    srand((unsigned int)time(NULL));  //seeding the randomize function

    cout<<"Enter n:";
    cin>>n;
    int m=(int)ceil(sqrt(n));  //definition of m
    indexnode** hash=inithash(m);  //hash initialized
    
    arrelem* arr=(arrelem*)malloc(n*sizeof(arrelem));     //point array formed
    
    initializepoints(arr,hash,m,n);  //allocate points
    
    rajasekaran(arr,hash,n,m);  //implementation of rajasekaran
}