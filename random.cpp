#include <iostream>
#include <stdlib.h>
#include <iomanip>

using namespace std;


int main(void)
{
    srand(time(NULL));
    int p = (rand()%8) + 2;
    int t = (rand()%10) + 2;
    int k = (rand()%10) + 2;
    int cof = (rand()%10);

    int n = p*k*t;

    double arr[n][n];

    int col = k* p;
    
    freopen ("IO/inputfile10.txt","w",stdout);

    printf("1\n");
    printf("%d\n",k);
    printf("%d\n",t);
    printf("%d\n",p);
    printf("%d\n",cof);

    for(int i =0; i< n;i++)
    {
        for(int j=0; j<=i;j++)
        {
            double random1 = (double)(rand()% 100)/99.0;
            //printf("%f.1 ", random1);

            if(i == j) arr[i][j] = 0.00;
            else 
            {
                arr[i][j] = random1;
                arr[j][i] = random1;
            }

            //cout << setprecision(2) << fixed << random1;

        }
        //printf("\n");
    }

    for(int i =0; i< n;i++)
    {
        for(int j=0; j< (n-1);j++)
        {
            //double random1 = (double)(rand()% 100)/99.0;
            //printf("%f.1 ", random1);

            // if(i == j) arr[i][j] = 0.00;
            // else 
            // {
            //     arr[i][j] = random1;
            //     arr[j][i] = random1;
            // }

            cout << setprecision(2) << fixed << arr[i][j]<<" ";

        }
        cout << setprecision(2) << fixed << arr[i][n-1];
        printf("\n");
    }

    /* code */
    return 0;
}
