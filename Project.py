import csv
import numpy as np
import pandas as pd
import ast
from time import time

def csvFileToAOA(fileName):
   with open(fileName) as csv_file:
      csv_reader = csv.reader(csv_file, delimiter=',')
      return list(csv_reader)

aoa_all = csvFileToAOA("SNP.example")

def mergepairs(aoa):
   new_aoa = []
   for row in aoa[1:len(aoa)]:
      new_row=[]

      for i in range(3,len(row),2):
         new_row.append(row[i] + row[i+1])
      new_aoa.append(new_row)
   return new_aoa

aoa = mergepairs(aoa_all)

def num_of_app(aoa):
   new_name = ["AA", "GG", "CC", "TT", "AG", "AC", "AT", "GT", "GT", "CT", ""]
   new_aoa = []
   new_aoa.append(new_name)

   for row in aoa:
      maxi = 0
      flag = 1
      new_row = [0,0,0,0,0,0,0,0,0,0,0]


      for i in range(len(row)):
         if(row[i] == "AA"):
            new_row[0] += 1
         if (row[i] == "GG"):
            new_row[1] += 1
         if (row[i] == "CC"):
            new_row[2] += 1
         if (row[i] == "TT"):
            new_row[3] += 1
         if (row[i] == "AG" or row[i] == "GA" ):
            new_row[4] += 1
         if (row[i] == "AC" or row[i] == "CA" ):
            new_row[5] += 1
         if (row[i] == "AT" or row[i] == "TA"):
            new_row[6] += 1
         if (row[i] == "GC" or row[i] == "CG"):
            new_row[7] += 1
         if (row[i] == "GT" or row[i] == "TG"):
            new_row[8] += 1
         if (row[i] == "CT" or row[i] == "TC" ):
            new_row[9] += 1
         if(row[i] == "00"):
            new_row[10] +=1
      for i in new_row:
         if(i>=maxi):
            maxi=i


      for i in range(len(new_row)):
         if (new_row[i] == maxi):
            new_row[i]+= new_row[10]
            new_row[10] = i

      new_aoa.append(new_row)
   return new_aoa

#print(num_of_app(aoa))

def new_matrix(aoa1):
   aoa = num_of_app(aoa1)
   new_matrix = []
   j = 0
   for row in aoa[1:]:
      A = ""
      a = ""
      index_maxi = ""
      flag =0
      for i in range(4):
         if (row[i] != 0):
            flag=1
            A = aoa[0][i]
            #print(A)
            break
      if(flag == 1):
         for i in range(4, 10):
            if (row[i] != 0):
               a = aoa[0][i]

               index_maxi = aoa[0][row[10]]

               break
      for row1 in aoa1[j:j+1]:
         new_row = []
         for i in range(len(row1)):
            if(row1[i] == "00"):
               row1[i] = index_maxi
            if(row1[i] == A):
               new_row.append(2)
               continue
            if(row1[i] == a or row1[i] == a[::-1]):
               new_row.append(1)
               continue
            else:
               new_row.append(0)
               continue
         new_matrix.append(new_row)
      j+=1
   #print(new_matrix)
   return new_matrix

def Equals(cow1,cow2):
   counter = 0
   for i in range(len(cow1)):
      if (cow1[i] == cow2[i]):
         counter+= 2
      else:
         if((cow1[i] == 2 and cow2[i] == 1) or (cow1[i] == 1 and cow2[i] == 2)):
            counter+= 1
            continue
         if((cow1[i] == 1 and cow2[i] == 0) or (cow1[i] == 0 and cow2[i] == 1)):
            counter+= 1
            continue

   return counter

def IBS(aoa_num):
   rowsSize = len(aoa[0])
   colSize = len(aoa[0])
   x = np.ones((rowsSize, colSize))

   f = open("IBS.csv", "w+")
   #m - (len(aoa_num[0]))
   for i in range((len(aoa_num[0]))):
      new_row = []

      X = [row[i] for row in aoa_num]
      #m (len(aoa_num[0]))
      for j in range(i,len(aoa_num[0])):
         Y = [row[j] for row in aoa_num]

         x[i][j] = Equals(X,Y)/(2*float(len(aoa_num)))
         x[j][i] = x[i][j]
         #new_row.append(Equals(X,Y)/(2*float(len(aoa_num))))
      #new_matrix.append(new_row)

   print(x)
   for e in x:
      for value in e:
         f.write(str(value))
         f.write(' ')
      f.write('\n')

   return x




def pc(aoa):
   new_p = []
   for row in aoa:
      counter = 0
      for i in range(len(row)):
         counter += row[i]/2

      new = counter / float((100))
      new_p.append(new)
   return new_p




#print(pc(new_matrix(aoa)))


def sumP(arr):
   sum=0
   for i in arr:
      sum += 2*i*(1-i)
   return sum
#print(sumP(pc(new_matrix(aoa))))

def IbdForPair(individualONE,individualTWO):

   counter = 0
   for i in range(len(individualONE)):
      counter+= (individualTWO[i] - individualONE[i])**2

   return counter

#print(num_of_app(aoa))

def IBD(aoa_num):
   rowsSize = len(aoa[0])
   colSize = len(aoa[0])
   x = np.ones((rowsSize, colSize))
   p = pc(aoa_num)
   counterP = sumP(p)



   f = open("IBD.csv", "w+")
   #len(aoa[0])
   for i in range(len(aoa_num[0])):
      new_row = []
      X = [row[i] for row in aoa_num]
      #len(aoa[0])
      for j in range(i,len(aoa_num[0])):
         Y = [row[j] for row in aoa_num]

         x[i][j] = 0.5 - (IbdForPair(X, Y))/(4*float(counterP))
         x[j][i] = x[i][j]
         #new_row.append(0.5 - (IbdForPair(X, Y))/(4*float(counterP)))
      #new_matrix.append(new_row)

   #print (new_matrix)
   print(x)
   for e in x:
      for value in e:
         f.write(str(value))
         f.write(' ')
      f.write('\n')

   return x



IBS(new_matrix(aoa))
#IBD(new_matrix(aoa))


def merge(aoa):
   new_aoa = []
   for row in aoa:
      new_row = []
      new_row.append(row[1])
      new_row.append(row[2])
      # len(row) - midelle
      for i in range(3, len(row), 2):
         new_row.append(row[i] + row[i + 1])
      new_aoa.append(new_row)
   return new_aoa

aoa_for_roh = merge(aoa_all)


def roh_individual(aoa):
   new_matrix = []
   f = open("ROH.csv", "w+")
   new_row = []
   new_row.append([row[0] for row in aoa])
   new_row.append([row[1] for row in aoa])
   for i in range(2, len(aoa[0])):
      new_row.append([row[i] for row in aoa])
   for i in range(2, len(new_row)):
      flagOne = 0
      new_cow = [0, 0, 0, 0, 0, 0]
      counter = 0
      denominator = 0
      if (new_row[i][0] == "AA" or new_row[i][0] == "CC" or new_row[i][0] == "TT" or new_row[i][0] == "GG"):
         flagOne = -1
      for j in range(1, len(new_row[0])):
         if (new_row[0][j] == new_row[0][flagOne + 1]):
            if (new_row[i][j] != "AA" and new_row[i][j] != "CC" and new_row[i][j] != "TT" and new_row[i][j] != "GG"):
               if (j != flagOne + 1):
                  x = (int(new_row[1][j - 1]) - int(new_row[1][flagOne + 1])) / float(1000000)
                  if (x < 4.0 and x > 2.0):
                     new_cow[0] += 1
                     counter += 3 * 1000
                     denominator += 1
                  if (x < 6.0 and x > 4.0):
                     new_cow[1] += 1
                     counter += 5 * 1000
                     denominator += 1
                  if (x < 8.0 and x > 6.0):
                     new_cow[2] += 1
                     counter += 7 * 1000
                     denominator += 1

                  if (x < 10.0 and x > 8.0):
                     new_cow[3] += 1
                     counter += 9 * 1000
                     denominator += 1
                  if (x > 10):
                     new_cow[4] += 1
                     counter += 10 * 1000
                     denominator += 1

                  flagOne = j
               else:
                  flagOne = j
            else:
               if (j == len(new_row[0]) - 1):
                  x = (int(new_row[1][j]) - int(new_row[1][flagOne + 1])) / float(1000000)
                  if (x < 4.0 and x > 2.0):
                     new_cow[0] += 1
                     counter += 3 * 1000
                     denominator += 1
                  if (x < 6.0 and x > 4.0):
                     new_cow[1] += 1
                     counter += 5 * 1000
                     denominator += 1
                  if (x < 8.0 and x > 6.0):
                     new_cow[2] += 1
                     counter += 7 * 1000
                     denominator += 1

                  if (x < 10.0 and x > 8.0):
                     new_cow[3] += 1
                     counter += 9 * 1000
                     denominator += 1
                  if (x > 10):
                     new_cow[4] += 1
                     counter += 10 * 1000
                     denominator += 1
         else:
            if (new_row[i][j] == "AA" or new_row[i][j] == "CC" or new_row[i][j] == "TT" or new_row[i][j] == "GG"):
               x = (int(new_row[1][j - 1]) - int(new_row[1][flagOne + 1])) / float(1000000)
               if (x < 4.0 and x > 2.0):
                  new_cow[0] += 1
                  counter += 3 * 1000
                  denominator += 1
               if (x < 6.0 and x > 4.0):
                  new_cow[1] += 1
                  counter += 5 * 1000
                  denominator += 1
               if (x < 8.0 and x > 6.0):
                  new_cow[2] += 1
                  counter += 7 * 1000
                  denominator += 1

               if (x < 10.0 and x > 8.0):
                  new_cow[3] += 1
                  counter += 9 * 1000
                  denominator += 1
               if (x > 10):
                  new_cow[4] += 1
                  counter += 10 * 1000
                  denominator += 1
            flagOne = j - 1

      new_cow[5] = counter / denominator
      new_matrix.append(new_cow)
   f.write(str(new_matrix))
   print(new_matrix)



roh_individual(aoa_for_roh)



def dis(cow1,cow2):
   distans = []
   for i in range(len(cow1)):
      distans.append(abs(int(cow2[i]) - int(cow1[i])))
   return distans


def distans(aoa,n):
   p = pc(aoa)

   rowsSize = len(aoa[0])
   colSize = len(aoa[0])
   x = np.ones((rowsSize, colSize))

   new_matrix = []
   f = open("DISTANS" + str(n) + "P.csv", "w+")
   #len(aoa[0])
   for i in range(len(aoa[0])):
      new_row = []
      u = 0
      X = [row[i] for row in aoa]
      #(len(aoa[0]))
      for j in range(i,len(aoa[0])):
         Y = [row[j] for row in aoa]
         #print(X)
         #print(Y)
         distans = dis(X,Y)
         #print(distans)

         u = 0
         for k in range(len(distans)-n+1):
            sum = 0
            for g in range(k,k+n):
               sum+=distans[g]
            u+= (4*n**2 - sum**2)*4*p[k]*(1-p[k])

         #new_row.append(u/float(len(distans)*4*n**2))
         x[i][j] = u/float((len(distans)-n+1)*4*n**2)
         x[j][i] = x[i][j]
      #new_matrix.append(new_row)
   # f.write("[")
   for e in x:
      for value in e:
         f.write(str(value))
         f.write(' ')
      f.write('\n')


   # f.write("]")
   # print(x)
   #x.to_csv('file.csv')
   #f.write(str(x))

   return x







#print len(aoa)

#y = np.genfromtxt('IBD.csv')
#y_df = pd.DataFrame(y)
#for index, row in y_df.describe().iterrows():
#   print([row[i] for i in range(100)])


#print(y_df.describe())





start = time()
#distans(new_matrix(aoa),1)

total = time()-start
print(total)
#print(total)


def normal(data):
   transX = np.copy(data.T)
   for j in range(100):
      for i in range(100):
         if np.var(transX[j]) == 0:
            continue
         data[i][j] = (transX[j][i] - np.average(transX[j])) / np.sqrt(np.var(transX[j]))
   print(data)

#normal(y)

