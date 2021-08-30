import numpy as np
import sys

def check():
    while True:
        try:
            m = input('Enter number of unknowns:')
            x = int(m)
            if x >= 2 and x <= 4:
                break
            else:
                print("Number less than 2 or greater than 4, try again")
        except ValueError:
            print("Must be a number, try again")
    return x

def gauss(l):
    g=l.shape[0]
    f=np.zeros((g,f+1))
    sln=np.zeros(g)

    if g == 2:
        c_1=float(input('Enter C_1:\n'))
        c_2=float(input('Enter C_2:\n'))

        f[0,:] =[l[0,:],c_1]
        f[1,:]=[l[1,:],c_2]

        for i in range(g):
            if f[i][i] == 0.0:
                sys.exit('Divide by zero detected!')

            for j in range(i+1, g):
                ratio = f[j][i]/f[i][i]

                for k in range(g+1):
                    f[j][k] = f[j][k] - ratio * f[i][k]

        sln[g-1] = f[g-1][g]/f[g-1][g-1]


        for i in range(g-2,-1,-1):
            sln[i] = f[i][g]

            for j in range(i+1,g):
                sln[i] = sln[i] - f[i][j]*sln[j]

            sln[i] = sln[i]/f[i][i]

        print('\nLinear system values: ')
        for i in range(g):
            print('X%d = %0.2f' %(i,sln[i]), end = '\t')

    if g == 3:
         d_1=float(input('Enter d_1:\n'))
         d_2=float(input('Enter d_2:\n'))
         d_3=float(input('Enter d_3:\n'))

         f[0,:]=[l[0,:],d_1]
         f[1,:]=[l[1,:],d_2]
         f[2,:]=[l[2,:],d_3]

         for i in range(g):
             if f[i][i] == 0.0:
                 sys.exit('Divide by zero detected!')

             for j in range(i+1, g):
                 ratio = f[j][i]/f[i][i]

                 for k in range(g+1):
                     f[j][k] = f[j][k] - ratio * f[i][k]

         sln[g-1] = f[g-1][g]/f[g-1][g-1]


         for i in range(g-2,-1,-1):
             sln[i] = f[i][g]

             for j in range(i+1,g):
                 sln[i] = sln[i] - f[i][j]*sln[j]

             sln[i] = sln[i]/f[i][i]

         print('\nLinear system values: ')
         for i in range(g):
             print('X%d = %0.2f' %(i,sln[i]), end = '\t')


    if g == 4:
             d_1=float(input('Enter d_1:\n'))
             d_2=float(input('Enter d_2:\n'))
             d_3=float(input('Enter d_3:\n'))
             d_4=float(input('Enter d_4:\n'))

             f[0,:]=[l[0,:],d_1]
             f[1,:]=[l[1,:],d_2]
             f[2,:]=[l[2,:],d_3]
             f[3,:]=[l[3,:],d_4]

             for i in range(g):
                 if f[i][i] == 0.0:
                     sys.exit('Divide by zero detected!')

                 for j in range(i+1, g):
                     ratio = f[j][i]/f[i][i]

                     for k in range(g+1):
                         f[j][k] = f[j][k] - ratio * f[i][k]

             sln[g-1] = f[g-1][g]/f[g-1][g-1]


             for i in range(g-2,-1,-1):
                 sln[i] = f[i][g]

                 for j in range(i+1,g):
                     sln[i] = sln[i] - f[i][j]*sln[j]

                 sln[i] = sln[i]/f[i][i]

             print('\nLinear system values: ')
             for i in range(g):
                 print('X%d = %0.2f' %(i,sln[i]), end = '\t')


#If the main determinant is zero the system of linear equations is either inconsistent or has infinitely many solutions try gauss
def cramer_2x2(m):

    d=np.linalg.det(m)
    if d == 0:
        gauss(m)

    else:
        c_1=float(input('Enter C_1:\n'))
        c_2=float(input('Enter C_2:\n'))

        det_x=np.linalg.det([[c_1,m[0,1]],[c_2,m[1,1]]])
        det_y=np.linalg.det([[m[0,0],c_1],[m[1,0],c_2]])

        sys_solution=list()
        sys_solution.append(det_x/d)
        sys_solution.append(det_y/d)
        print("Linear system values:",*sys_solution)

def cramer_3x3(p):
    d=np.linalg.det(p)
    if d == 0:
        gauss(p)
    else:
        d_1=float(input('Enter d_1:\n'))
        d_2=float(input('Enter d_2:\n'))
        d_3=float(input('Enter d_3:\n'))

        sys_val=[d_1,d_2,d_3]
        x=p.copy()
        x[:,0]=sys_val
        det_x=np.linalg.det(x)

        y=p.copy()
        y[:,1]=sys_val
        det_y=np.linalg.det(y)

        z=p.copy()
        z[:,2]=sys_val
        det_z=np.linalg.det(z)

        sys_solution=list()
        sys_solution.append(det_x/d)
        sys_solution.append(det_y/d)
        sys_solution.append(det_z/d)
        print("Linear system values:",*sys_solution)


def cramer_4x4(j):
    d=np.linalg.det(j)
    if d == 0:
        gauss(j)
    else:
        d_1=float(input('Enter d_1:\n'))
        d_2=float(input('Enter d_2:\n'))
        d_3=float(input('Enter d_3:\n'))
        d_4=float(input('Enter d_4:\n'))

        sys_val=[d_1,d_2,d_3,d_4]
        x=j.copy()
        x[:,0]=sys_val
        det_x=np.linalg.det(x)

        y=j.copy()
        y[:,1]=sys_val
        det_y=np.linalg.det(y)

        z=j.copy()
        z[:,2]=sys_val
        det_z=np.linalg.det(z)

        w=j.copy()
        w[:,3]=sys_val
        det_w=np.linalg.det(w)

        sys_solution=list()
        sys_solution.append(det_x/d)
        sys_solution.append(det_y/d)
        sys_solution.append(det_z/d)
        sys_solution.append(det_w/d)
        print("Linear system values:",*sys_solution)


#number of unknowns
n=check()
#coefficients of unknowns
matrix_of_coefficients=np.zeros((n,n))

print('Enter the coefficients of unknowns:\n')
for i in range(n):
    for j in range(n):
        matrix_of_coefficients[i][j] = float(input( 'coefficients['+str(i)+']['+ str(j)+']='))
        while matrix_of_coefficients[i][j] == 0:
            print('Coefficients cannot be zero.\n')
            matrix_of_coefficients[i][j] = float(input( 'coefficients['+str(i)+']['+ str(j)+']='))


print("This is the matrix of coefficients:\n",matrix_of_coefficients,"\nDo you want to change any element?\n")
response=['Yes','No']
r=input("If so input Yes else No.")

while r not in response:
    print("Unrecognized input!")
    r=input("If so input Yes else No.")

row=0
col=0

if r == response[0]:
    row=int(input('Enter row number:\n'))
    col=int(input('Column number:\n'))

    while  (row > n or row < 0):
          print('Row number cannot be',row)
          row=int(input('Enter row number:\n'))
    while (col > n or col < 0):
          print('Column number cannot be',col)
          col=int(input('Column number:\n'))


    matrix_of_coefficients[row][col]=float(input('Enter new value:\n'))

    print('Changed matrix:\n',matrix_of_coefficients)

if matrix_of_coefficients.shape[0] == 2:
    cramer_2x2(matrix_of_coefficients)

elif matrix_of_coefficients.shape[0] == 3:
    cramer_3x3(matrix_of_coefficients)

else:
    cramer_4x4(matrix_of_coefficients)
