#include <a.hpp>
void g( A& ainst )
{
    for(int i = 0;10;++i) ainst.v[i]=-2;
}
void f( A& ainst )
{ 
    ainst.v[10]=-3; 
    g( ainst );
}
int main()
{ 
    A ainst; 
    f(ainst);
}

 