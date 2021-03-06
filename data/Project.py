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


#IBD(new_matrix(aoa))


def merge(aoa):
   new_aoa = []
   for row in aoa[1:len(aoa)]:
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



#roh_individual(aoa_for_roh)



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

def sum_p(individualONE, individualTWO):
    counter = 0
    for i in range(len(individualONE)):
        if(individualONE[i][0] == individualONE[i][1] and individualTWO[i][0] == individualTWO[i][1] and individualTWO[i][0] == individualONE[i][0]):
            counter+=1
        elif(individualONE[i] != individualTWO[i] and individualONE[i][0] == individualONE[i][1] and individualTWO[i][1] == individualTWO[i][0]):
            counter+=0
        else:
            counter+=0.5

    return counter

def P_Homozygosity(aoa):
   rowsSize = len(aoa[0])
   colSize = len(aoa[0])
   x = np.ones((rowsSize, colSize))
   f = open("P_HOMO.csv", "w+")
   #len(aoa[0])
   for i in range(len(aoa[0])):
      new_row = []
      cow1 = [row[i] for row in aoa]
      #len(aoa[0])
      for j in range(i,len(aoa[0])):
         cow2 = [row[j] for row in aoa]
         x[i][j] = sum_p(cow1, cow2)/float(len(cow1))
         x[j][i] = x[i][j]

   for e in x:
      for value in e:
         f.write(str(value))
         f.write(' ')
      f.write('\n')

   return x

#P_Homozygosity(aoa)

def ArrayEquals(cow1,cow2):
   distan = np.zeros((len(cow1)))

   for i in range(len(cow1)):
      if (cow1[i] == cow2[i]):
         distan[i] = 2
         continue
      else:
         if((cow1[i] == 2 and cow2[i] == 0) or (cow1[i] == 0 and cow2[i] == 2)):
            distan[i] = 0
            continue
         else:
            distan[i] = 1
            continue
   print(distan)
   return distan

def runs_ibd2_ibd1(distans, matrix, numcow1, numcow2):
   cromozom = matrix[0]
   print(cromozom)
   snp_cow1 = matrix[numcow1]
   #print(snp_cow1)
   snp_cow2 = matrix[numcow2]
   #print(snp_cow2)
   i = 0
   start = 0
   end = 0
   sum = 0
   new_row = [0,0,0,0]
   while (i<len(distans)):

      if (distans[i] == 0 and cromozom[end+1] == cromozom[end] and i+1<len(distans)):
         start = i
         end = i
         while (i+1 < len(distans) and cromozom[end+1] == cromozom[end]):
            if (distans[i] != 0):
               if ((i + 1 < len(distans) and distans[i] == 2 and distans[i + 1] == 0 and cromozom[i+1] == cromozom[end])):
                  distans[i] = 0
                  i = i + 1
                  end = i

               else:
                  break
            else:
               if((i+2<len(distans) and distans[i+1] != 2 and distans[i+2] !=0) or (i+2 == len(distans) and distans[i+1] != 0 and start == i)) :
                  break
               else:
                  i = i + 1
                  end = i


      i = i + 1
      start = 0
      end = 0
   print("dis after 0:")
   print(distans)
   i = 0

   while (i < len(distans)):
      if (distans[i] != 2 and cromozom[end + 1] == cromozom[end] and i + 1 < len(distans)):
         start = i
         end = i
         while (i + 1 < len(distans) and cromozom[end + 1] == cromozom[end]):
            if (distans[i] == 2):
               if ((i + 1 < len(distans) and distans[i] == 2 and distans[i + 1] != 2 and cromozom[i + 1] == cromozom[end])):
                  distans[i] = 1
                  i = i + 1
                  end = i

               else:
                  break
            else:
               if ((i + 2 < len(distans) and distans[i + 1] != 2 and distans[i + 2] == 2) or (
                       i + 2 == len(distans) and distans[i + 1] == 2 and start == i)):
                  break
               else:
                  i = i + 1
                  end = i

      i = i + 1
      start = 0
      end = 0
   print("dis after 121:")
   print(distans)
   i = 0

   while (i<len(distans)):
      S1 = ""
      S2 = ""
      if (distans[i] == 2  and i+1<len(distans)  and cromozom[end+1] == cromozom[end]):
         start = i
         end = i

         while (i+1 < len(distans) and cromozom[i+1] == cromozom[end]):

            if (distans[i] != 2):
               if ((i + 1 < len(distans) and distans[i] != 2 and distans[i + 1] == 2 and cromozom[i+1] == cromozom[end])  ):
                  print("kkk")
                  distans[i] = 2
                  snp_cow1[i] = "AA"
                  snp_cow2[i] = "AA"
                  i = i + 1
                  end = i
               else:
                  break
            else:
               if((i+2 < len(distans) and distans[i+1] != 2 and distans[i+2] != 2  and cromozom[i+1] != cromozom[end]) or
                  (i+2 < len(distans) and distans[i + 1] != 2 and distans[i + 2] != 2 and  cromozom[i+1] == cromozom[end]) or

                  ( i+2 == len(distans) and distans[i+1] != 2 and start == i)) :
                  distans[start]=1
                  break
               else:
                  i = i + 1
                  end = i


      elif((i+1<len(distans) and distans[i] == 2 and cromozom[end+1] != cromozom[end]) or (distans[i] == 2 and i+1 == len(distans))):
         distans[i] =1




      if(end != start):
         print("start:" + str(start))
         print("end:" +str(end))
         for k in range(start, end+1):
            S1 += str(snp_cow1[k][0])
            S2 += str(snp_cow1[k][1])
         print(S1)
         print(S2)
         if (S1 != ""):
            if (S1 == S2):
               counter = int(matrix[1][end]) - int(matrix[1][start])
               # print("conter:::")
               # print(counter)
               # 1000000
               if (counter < 50000):
                  for k in range(start, end + 1):
                     distans[k] = 1
               else:
                  for k in range(start, end + 1):
                     distans[k] = -1
                  new_row[0] += 1
                  sum += counter
            else:
               # print("end crom")
               # print(int(matrix[1][end]))
               # print("start crom")
               # print(matrix[1][start])
               counter = 0.5 * (int(matrix[1][end]) - int(matrix[1][start]))
               # print("conter:")
               # print(counter)
               # 500000
               if (counter < 50000):
                  for k in range(start, end + 1):
                     distans[k] = 1
               else:
                  for k in range(start, end + 1):
                     distans[k] = -1
                  new_row[1] += 1
                  sum += counter
      elif(i==start and start != 0):
         distans[start] = 1

      i = i + 1
      start = 0
      end = 0
   print("sum after 2:")
   print(sum)
   print("dis after 2:")
   print(distans)
   i = 0


   while (i<len(distans)):

      if (distans[i] == 0 and cromozom[end+1] == cromozom[end] and i+1<len(distans)):
         start = i
         end = i
         while (i+1 < len(distans) and cromozom[end+1] == cromozom[end]):
            if (distans[i] != 0):
               if ((i + 1 < len(distans) and distans[i] == 1 and distans[i + 1] == 0 and cromozom[i+1] == cromozom[end])):
                  distans[i] = 0
                  i = i + 1
                  end = i

               else:
                  break
            else:
               if((i+2<len(distans) and distans[i+1] != 1 and distans[i+2] !=0) or (i+2 == len(distans) and distans[i+1] != 0 and start == i)) :
                  break
               else:
                  i = i + 1
                  end = i
      elif((i+1<len(distans) and distans[i] == 2 and cromozom[end+1] != cromozom[end]) or (distans[i] == 2 and i+1 == len(distans))):
         distans[i] =1

      i = i + 1
      start = 0
      end = 0
   print("sum after 0:")
   print(sum)
   print("dis after 0:")
   print(distans)
   i = 0







   while (i<len(distans)):
      S1 = ""
      S2 = ""
      S3 = ""


      if ((distans[i] == 1 and i+1<len(distans) and cromozom[end+1] == cromozom[end])):
         start = i
         end = i
         while (i+1 < len(distans) and cromozom[end+1] == cromozom[end]):
            if (distans[i] != 1):
               if ((i + 1 < len(distans) and distans[i] == 0 and distans[i + 1] == 1 and cromozom[i+1] == cromozom[end])):
                  distans[i] = 1
                  snp_cow1[i] = "AA"
                  snp_cow2[i] = "AA"
                  i = i + 1
                  end = i
               else:
                  break
            else:

               if((i+2<len(distans) and distans[i+1] != 1 and distans[i+2] !=1 and cromozom[i+1] != cromozom[end] )
                       or (i+2<len(distans) and distans[i+1] != 1 and distans[i+2] !=1 and cromozom[i+1] == cromozom[end]) or
                       (i+2 == len(distans) and distans[i+1] != 1 and start == i)) :
                  distans[start]=1
                  break
               else:

                  i = i + 1
                  end = i

      elif((i+1<len(distans) and distans[i] == 1 and cromozom[end+1] != cromozom[end]) or (distans[i] == 1 and i+1 == len(distans))):
         distans[i] =0


      if(start != end):
         print("start:" + str(start))
         print("end:" + str(end))
         for k in range(start, end + 1):
            if (snp_cow1[k][0] == snp_cow2[k][0]):
               S1 += str(snp_cow1[k][0])
               S2 += str(snp_cow1[k][1])
               S3 += str(snp_cow2[k][1])
            elif (snp_cow1[k][0] == snp_cow2[k][1]):
               S1 += str(snp_cow1[k][0])
               S2 += str(snp_cow1[k][1])
               S3 += str(snp_cow2[k][0])
            elif (snp_cow1[k][1] == snp_cow2[k][0]):
               S1 += str(snp_cow1[k][1])
               S2 += str(snp_cow1[k][0])
               S3 += str(snp_cow2[k][1])
            elif (snp_cow1[k][1] == snp_cow2[k][1]):
               S1 += str(snp_cow1[k][1])
               S2 += str(snp_cow1[k][0])
               S3 += str(snp_cow2[k][0])
         # print(S1)
         # print(S2)
         # print(S3)
         if (S1 != ""):
            if (S1 == S2 or S1 == S3):
               # print("end crom:")
               # print((matrix[1][end]))
               # print("start crom:")
               # print((matrix[1][start]))
               counter = 0.5 * (int(matrix[1][end]) - int(matrix[1][start]))
               # print("the counter:")
               # print(counter)
               # 500000
               if (counter < 50000):
                  for k in range(start, end + 1):
                     distans[k] = 0
               else:
                  for k in range(start, end + 1):
                     distans[k] = -1
                  new_row[2] += 1
                  sum += counter
            else:
               # print("end crom:")
               # print((matrix[1][end]))
               # print("start crom:::")
               # print((matrix[1][start]))
               counter = 0.25 * (int(matrix[1][end]) - int(matrix[1][start]))
               # print("the counter:")
               # print(counter)
               # 250000
               if (counter < 50000):
                  for k in range(start, end + 1):
                     distans[k] = 0
               else:
                  for k in range(start, end + 1):
                     distans[k] = -1
                  new_row[3] += 1
                  sum += counter
      elif(start==i and distans[start]!=-1):
         distans[start] = 0



      i = i + 1
      start = 0
      end = 0

   print("sum:")
   print(sum)
   print("dis :")
   print(distans)
   print(new_row)
   return sum






def matrix_Equals(aoa, aoa_all):
   rowsSize = len(aoa[0])
   colSize = len(aoa[0])
   x = np.ones((rowsSize, colSize))
   f = open("run_ibd2and1.csv", "w+")
   new_row = []
   new_row.append([row[0] for row in aoa_all])
   new_row.append([row[1] for row in aoa_all])
   # len(aoa[0])
   for i in range(2, len(aoa_all[0])):
      new_row.append([row[i] for row in aoa_all])

   # len(aoa[0])
   for i in range(1):
      cow1 = [row[i] for row in aoa]
      # len(aoa[0])
      for j in range(55,56):
         cow2 = [row[j] for row in aoa]
         dis = ArrayEquals(cow1, cow2)
         #dis = [2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ]
         x[i][j] = runs_ibd2_ibd1(dis, new_row, i + 2, j + 2)
         x[j][i] = x[i][j]
         # new_row.append(0.5 - (IbdForPair(X, Y))/(4*float(counterP)))
      # new_matrix.append(new_row)

   return x

   # print (new_matrix)
   # print(x)


matrix_Equals(new_matrix(aoa), aoa_for_roh)

