#include <stdio.h>
#include <chrono>
#include <random>

using namespace std;

int main()
{
    /*for (int i=0;i<5;i++)
    {
        int a=(i-1+5)%5;
        int b=(i+1)%5;
        printf("a: %i\n", a);
        printf("b: %i\n", b);
    }*/
    //printf("%i", (4+1+5)%5);
    /*unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
    mt19937_64 generator(seed1);   
    uniform_int_distribution<int> i_distribution(-1,1);
    for (int i=0;i<10;i++)
    {
        int a=i_distribution(generator);
        printf("%i\n", a);
    }*/

    for (int k=0;k<=1000;k++)
    {
        if (k%10==0)
        {
            printf("%i\n", k);
        }
    }

    return 0;
}