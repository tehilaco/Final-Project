import urllib.request
import csv
import numpy as np
import pandas as pd

def csvFileToAOA(fileName):
   with open(fileName) as csv_file:
      csv_reader = csv.reader(csv_file, delimiter=',')
      return list(csv_reader)

aoa = csvFileToAOA("moran_snp.csv")

def urlcromANDlocatin():
    url = "https://www.animalgenome.org/repository/cattle/Illinoi_Beever_Project.2012/SNP_Map_13oct2009.txt"
    file = urllib.request.urlopen(url)
    f = open("CromLocation.csv", "w+")
    for line in file:
	    decoded_line = line.decode("utf-8")
	    f.write(decoded_line)
    with open("CromLocation.csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        return list(csv_reader)

CromLocation = urlcromANDlocatin()

#print(CromLocation)

def Crom_location(aoa):
   new_aoa = []
   for row in aoa[2:len(aoa):2]:
      new_row=[]
      for i in range(1,2):
         new_row.append(row[i])
      for i in range(2, 4):
         if(row[i] != 'X'):
            new_row.append(int(row[i]))
         else:
            row[i] = 30
            new_row.append(int(row[i]))
      new_aoa.append(new_row)
   return new_aoa

CromAndlocation = Crom_location(CromLocation)
sorter = lambda x: (x[1], x[2])
new_list = sorted(CromAndlocation, key= sorter)

def merge(aoa , name, num_crom , location):
   new_row = []
   for row in aoa[1:len(aoa)]:
      if(row[0] == name):
         new_row.append(name)
         new_row.append(num_crom)
         new_row.append(location)
         for i in range(1, len(row)):
            new_row.append(row[i])
         return new_row
      else:
         continue
   return new_row

def merge_name(aoa,sortbycrom):
   new_aoa = []
   new_row = []
   new_row.append('MARKER_NAME')
   new_row.append('Chr')
   new_row.append('MapInfo')
   for row in aoa[0:1]:
         for i in range(1,len(row)):
            new_row.append(row[i])
   new_aoa.append(new_row)
   for index in sortbycrom[0:len(sortbycrom)]:
      x = merge(aoa,index[0],index[1],index[2])
      if(x != []):
         new_aoa.append(x)
   return new_aoa

aoa_all = merge_name(aoa,new_list)
#print(aoa_all)

def mergepairs(aoa):
   new_aoa = []
   for row in aoa[1:len(aoa)]:
      new_row=[]
      for i in range(3,len(row),2):
         new_row.append(row[i] + row[i+1])
      new_aoa.append(new_row)
   return new_aoa

aoa = mergepairs(aoa_all)
#print(aoa)

def num_of_app(aoa):

   new_name = ["AA", "GG", "CC", "TT", "AG", "AC", "AT", "GC", "GT", "CT", ""]
   new_aoa = []
   new_aoa.append(new_name)
   for row in aoa:
      maxi = 0
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
         if(row[i] == '  ' or row[i] == ' '):
            new_row[10] +=1


      for i in new_row[0:9]:
         if(i>=maxi):
            maxi=i

      for i in range(len(new_row)):
         if (new_row[i] == maxi):
            new_row[i]+= new_row[10]
            new_row[10] = i

      new_aoa.append(new_row)
   return new_aoa

def select_pair(aoa,Cow,Sire):
   cow = -2
   sire = -2
   counter = -2

   for row in aoa[0]:

      if row != "":
         counter += 1
      if row == str(Cow):
         cow = counter
      elif row == str(Sire):
         sire = counter

   return cow,sire



def Equals(cow1,cow2):
   counter = 0
   for i in range(len(cow1)):
      if (cow1[i] == cow2[i]):
         counter += 2
      else:
         if((cow1[i] == 2 and cow2[i] == 1) or (cow1[i] == 1 and cow2[i] == 2)):
            counter += 1
            continue
         if((cow1[i] == 1 and cow2[i] == 0) or (cow1[i] == 0 and cow2[i] == 1)):
            counter += 1
            continue

   return counter

def IBS(aoa_num , index ):
   i = index[0]
   j = index[1]
   cow = [row[i-2] for row in aoa_num]
   sire = [row[j-2] for row in aoa_num]

   x = Equals(cow,sire)/(2*float(len(aoa_num)))

   return x

def sum_p(cow1,cow2):
   counter = 0
   for i in range(len(cow1)):
      if (cow1[i] == cow2[i] and cow1[i] == 2):
         counter += 1
      elif(cow1[i] == cow2[i] and cow1[i] == 0):
         counter += 1
      elif((cow1[i] == 2 and cow2[i] == 0) or (cow1[i] == 0 and cow2[i] == 2)):
            counter += 0
      else:
         counter += 0.5

   return counter

def P_Homozygosity(aoa , index ):
   i = index[0]
   j = index[1]
   cow1 = [row[i-2] for row in aoa]
   cow2 = [row[j-2] for row in aoa]
   x = sum_p(cow1, cow2) / float(len(cow1))
   return x

def pc(aoa):
   new_p = []
   for row in aoa:

      counter = 0
      for i in range(len(row)):
         counter += row[i]/2

      new = counter / float(2286)
      new_p.append(new)
   #print(new_p)
   return new_p

def sumP(arr):
   sum=0
   for i in arr:
      sum += 2*i*(1-i)
   return sum

def IbdForPair(individualONE,individualTWO):
   counter = 0
   for i in range(len(individualONE)):
      counter+= (individualTWO[i] - individualONE[i])**2
   return counter

def IBD(aoa_num , index):
   i = index[0]
   j = index[1]
   p = pc(aoa_num)
   counterP = sumP(p)
   X = [row[i-2] for row in aoa_num]
   Y = [row[j-2] for row in aoa_num]

   x = 0.5 - (IbdForPair(X, Y))/(4*float(counterP))
   return x

def new_matrix(aoa1):

   aoa = num_of_app(aoa1)

   new_matrix = []
   j = 0
   for row in aoa[1:]:
      A = ""
      a = ""
      index_maxi = aoa[0][row[10]]
      flag =0
      for i in range(4):
         if (row[i] != 0):
            flag=1
            A = aoa[0][i]

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
            if(row1[i] == "  "):
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

   return new_matrix

#print(new_matrix(aoa))

def counter_ibs(cow1,cow2 , num):
   counter = 0
   if(num == 0):
      for i in range(len(cow1)):
         if ((cow1[i] == 2 and cow2[i] == 0) or (cow1[i] == 0 and cow2[i] == 2)):
            counter += 1
         else:
            counter += 0
   elif(num == 1):
      for i in range(len(cow1)):
         if ((cow1[i] == 2 and cow2[i] == 1) or (cow1[i] == 1 and cow2[i] == 2)
         or (cow1[i] == 1 and cow2[i] == 0) or (cow1[i] == 0 and cow2[i] == 1)):
            counter += 1
         else:
            counter += 0
   else:
      for i in range(len(cow1)):
         if (cow1[i] == cow2[i]):
            counter += 1
         else:
            counter += 0

   return counter

def pq00(aoa):
   counter = 0
   for row in aoa:
      p = 0
      q = 0
      for i in range(len(row)):
         p += row[i]
      q = 2*len(row) - p
      x = p / float(2 * len(row))
      y = q / float(2 * len(row))
      counter += 2*x*x*y*y
   return counter
def pq10(aoa):
   new_p = []
   new_q = []
   counter = 0
   for row in aoa:
      p = 0
      q = 0
      for i in range(len(row)):
         p += row[i]
      q = 2*len(row) - p
      x = p / float(2 * len(row))
      y = q / float(2 * len(row))
      counter += 4*x*x*x*y + 4*y*y*y*x
   return counter
def pq11(aoa):

   counter = 0
   for row in aoa:
      p = 0
      q = 0
      for i in range(len(row)):
         p += row[i]
      q = 2*len(row) - p
      x = p / float(2 * len(row))
      y = q / float(2 * len(row))
      counter += 2*x*x*y + 2*y*y*x
   return counter
def pq21(aoa):
   counter = 0
   for row in aoa:
      p = 0
      q = 0
      for i in range(len(row)):
         p += row[i]
      q = 2*len(row) - p
      x = p / float(2 * len(row))
      y = q / float(2 * len(row))
      counter +=  x*x*x + y*y*y + x*x*y + y*y*x
   return counter
def pq22(aoa):

   counter = 0
   for row in aoa:
      p = 0
      for i in range(len(row)):
         p += row[i]
      counter += 1
   return counter
def pq20(aoa):

   counter = 0
   for row in aoa:
      p = 0
      q = 0
      for i in range(len(row)):
         p += row[i]
      q = 2*len(row) - p
      x = p / float(2*len(row))
      y = q / float(2*len(row))
      counter += x*x*x*x + y*y*y*y + 4*x*x*y*y
   return counter


def P_IBS0(aoa_num , index ):
   i = index[0]
   j = index[1]
   cow = [row[i-2] for row in aoa_num]
   sire = [row[j-2] for row in aoa_num]
   N0 = counter_ibs(cow,sire , 0)
   return N0
def P_IBS1(aoa_num , index ):
   i = index[0]
   j = index[1]
   cow = [row[i-2] for row in aoa_num]
   sire = [row[j-2] for row in aoa_num]
   N1 = counter_ibs(cow,sire , 1)
   return N1
def P_IBS2(aoa_num , index ):
   i = index[0]
   j = index[1]
   cow = [row[i-2] for row in aoa_num]
   sire = [row[j-2] for row in aoa_num]

   N2 = counter_ibs(cow,sire , 2)

   return N2

#print(new_matrix(aoa))

#pANDq(new_matrix(aoa))


def new_ibd():
   df = pd.read_excel('Moran_cows_if.xlsx')
   new_ibs = []

   x = new_matrix(aoa)
   z = pq00(x)
   G = pq10(x)
   h = pq11(x)
   j = pq21(x)
   k = pq22(x)
   f = pq20(x)
   for i in range(len(df["Cow"])):
      y = select_pair(aoa_all, df["Dam"][i], df["Sire"][i])
      Z_0 = P_IBS0(x, y) / z
      Z_1 = (P_IBS1(x, y) - Z_0*G) / h
      Z_2 = (P_IBS2(x,y) - Z_0*f - Z_1*j) / k
      new_ibs.append(Z_2)

   df["IBD_new"] = new_ibs
   df.to_excel("new.xlsx")

#new_ibd()



def threeMadadim():
   x = new_matrix(aoa)
   df = pd.read_excel('Moran_cows_if.xlsx')
   new_ibs = []
   new_ibd = []
   new_ha = []
   i=0
   for i in range(len(df["Cow"])):
      print(i)
      y = select_pair(aoa_all, df["Dam"][i], df["Sire"][i])
      ibs_list = IBS(x, y)
      new_ibs.append(ibs_list)
      ha_list = P_Homozygosity(x, y)
      new_ha.append(ha_list)
      ibd_list = IBD(x, y)
      new_ibd.append(ibd_list)


   df["IBS"] = new_ibs
   df["IBD"] = new_ibd
   df["HA"] = new_ha

   df.to_excel("treeMA.xlsx")

#threeMadadim()

def location_dam(aoa,Cow):
   cow = -3
   counter = -3
   for row in aoa[0]:
      if row == str(Cow):
         cow = counter
      if row != "":
         counter += 1
         continue


   return cow


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
   #print(distan)
   return distan

def roh_individual(aoa,aoa_num, k):
   new_row = []
   new_row.append([row[0] for row in aoa])
   new_row.append([row[1] for row in aoa])
   new_row.append([row[k] for row in aoa_num])
   print(new_row[0])
   print(new_row[1])
   flagOne = 0
   counter = 0
   end = 0
   i=0
   print(new_row[2])
   
   while (i<len(new_row[2])):
      if (new_row[2][i] == 1 and new_row[0][end+1] == new_row[0][end] and i+1<len(new_row[2])):
         start = i
         end = i
         while (i+1 < len(new_row[2]) and new_row[0][end+1] == new_row[0][end]):
            if (new_row[2][i] != 1):
               if ((i + 1 < len(new_row[2]) and new_row[2][i + 1] == 1 and new_row[0][i+1] == new_row[0][end])):
                  new_row[2][i] = 1
                  i = i + 1
                  end = i
               else:
                  break
            else:
               if((i+2<len(new_row[2]) and new_row[2][i+1] == 1) or (i+2 == len(new_row[2]) and start == i)) :
                  break
               else:
                  i = i + 1
                  end = i
      i = i + 1

      end = 0
   i = 0
   
   
   while (i<len(new_row[2])):
      if (new_row[2][i] != 1  and i+1<len(new_row[2]) and new_row[0][end+1] == new_row[0][end]):

         start = i
         end = i
         while (i+1 < len(new_row[2]) and new_row[0][i+1] == new_row[0][end]):

            if (new_row[2][i] == 1):
               if ((i + 1 < len(new_row[2]) and new_row[2][i] == 1 and new_row[2][i + 1] != 1 and new_row[0][i+1] == new_row[0][end])):
                  new_row[2][i] = 2
                  i = i + 1
                  end = i
               else:
                  break
            else:

               if((i+2 < len(new_row[2]) and new_row[2][i+1] == 1 and new_row[2][i+2] == 1  and new_row[0][i+1] != new_row[0][end]) or
                  (i+2 < len(new_row[2]) and new_row[2][i + 1] == 1 and new_row[2][i + 2] == 1 and  new_row[0][i+1] == new_row[0][end]) or
                  (i+2 == len(new_row[2]) and new_row[2][i+1] == 1 and start == i)) :
                  break
               else:

                  i = i + 1
                  end = i
      if(i != 0 and new_row[2][i-1] != 1 and i+1<len(new_row[2]) and new_row[2][i] == 1 and new_row[0][i] == new_row[0][i-1] and new_row[0][i] != new_row[0][i+1]):
         new_row[2][i]  = 2
      i = i + 1
      end = 0
   print(new_row[2])

   if (new_row[2][0] == 2 or new_row[2][0] == 0):
      flagOne = -1
   for j in range(1, len(new_row[0])):
      if (new_row[0][j] == new_row[0][flagOne + 1]):
         if (j+1 < len(new_row[2]) and new_row[2][j] == 1):
            if (j != flagOne + 1):
               print(new_row[1][j - 1])
               print(new_row[1][flagOne + 1])
               x = (int(new_row[1][j - 1]) - int(new_row[1][flagOne + 1])) / float(1000000)
               if (x >= 1):
                   counter += x
               flagOne = j
            else:
               flagOne = j
         else:
            if (j == len(new_row[0]) - 1 or new_row[0][j] != new_row[0][j + 1]):
               print(new_row[1][j])
               print(new_row[1][flagOne + 1])
               x = (int(new_row[1][j]) - int(new_row[1][flagOne + 1])) / float(1000000)
               if ( x >= 1):
                  counter += x
               flagOne = j
      else:

         if (new_row[2][j] == 0 or new_row[2][j] == 2):
            print(new_row[1][j - 1])
            print(new_row[1][flagOne + 1])
            x = (int(new_row[1][j - 1]) - int(new_row[1][flagOne + 1])) / float(1000000)
            if (x >= 1):
                counter += x
         flagOne = j - 1

   return counter

def max_show(aoa):
   new_list = []
   x = num_of_app(aoa)
   for row in x[1:len(x)]:
      new_list.append(x[0][row[len(row)-1]])
   return new_list

def merge_max_gen(aoa,arr):
   new_aoa = []
   counter = 0
   for row in aoa[1:len(aoa)]:
      new_row = []
      new_row.append(row[1])
      new_row.append(row[2])
      # len(row) - midelle
      for i in range(3, len(row) , 2):
         if (row[i] != ' '):
            new_row.append(row[i] + row[i+1])
         else:
            new_row.append(arr[counter])
      new_aoa.append(new_row)
      counter += 1
   return new_aoa

aoa_for_roh = merge_max_gen(aoa_all,max_show(aoa))

def expect_roh(distans, matrix):
   cromozom = matrix[0]
   snp_cow1 = matrix[2]
   snp_cow2 = matrix[3]
   print(cromozom)
   i = 0
   start = 0
   end = 0
   sum = 0
   new_row = [0,0,0,0]
   print(distans)
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
   i = 0
   print(distans)
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
   i = 0
   print(distans)
   while (i<len(distans)):
      S1 = ""
      S2 = ""
      if (distans[i] == 2  and i+1<len(distans)  and cromozom[end+1] == cromozom[end]):
         start = i
         end = i
         while (i+1 < len(distans) and cromozom[i+1] == cromozom[end]):
            if (distans[i] != 2):
               if ((i + 1 < len(distans) and distans[i] != 2 and distans[i + 1] == 2 and cromozom[i+1] == cromozom[end])):
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
                  distans[i+1] = 2
                  snp_cow1[i+1] = "AA"
                  snp_cow2[i+1] = "AA"
                  i = i + 1
                  end = i
      elif((i+1<len(distans) and distans[i] == 2 and cromozom[end+1] != cromozom[end]) or (distans[i] == 2 and i+1 == len(distans))):
         distans[i] =1

      if(end != start):

         for k in range(start, end+1):
            S1 += str(snp_cow1[k][0])
            S2 += str(snp_cow1[k][1])

         if (S1 != ""):
            if (S1 == S2):
               counter = (int(matrix[1][end]) - int(matrix[1][start])) / float(1000000)
               if (counter < 1.0):
                  for k in range(start, end + 1):
                     distans[k] = 1
               else:
                  for k in range(start, end + 1):
                     distans[k] = -1
                  new_row[0] += 1
                  sum += counter
            else:
               counter = (int(matrix[1][end]) - int(matrix[1][start])) / float(1000000)
               if (counter < 1.0):
                  for k in range(start, end + 1):
                     distans[k] = 1
               else:
                  for k in range(start, end + 1):
                     distans[k] = -1
                  new_row[1] += 1
                  sum += (0.5 * counter)

      elif(i==start and start != 0):
         distans[start] = 1

      i = i + 1
      start = 0
      end = 0
   i = 0
   print(distans)
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
   i = 0
   print(distans)
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
         distans[i] = 0

      if(start != end):

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
         if (S1 != ""):
            if (S1 == S2 or S1 == S3):
               counter =  (int(matrix[1][end]) - int(matrix[1][start])) / float(1000000)
               if (counter <1.0):
                  for k in range(start, end + 1):
                     distans[k] = 0
               else:
                  for k in range(start, end + 1):
                     distans[k] = -1
                  new_row[2] += 1
                  sum += (0.5 * counter)
            else:
               counter =  (int(matrix[1][end]) - int(matrix[1][start])) / float(1000000)

               if (counter < 1.0):
                  for k in range(start, end + 1):
                     distans[k] = 0
               else:
                  for k in range(start, end + 1):
                     distans[k] = -1
                  new_row[3] += 1
                  sum += (counter * 0.25)

      elif(start==i and distans[start]!=-1):
         distans[start] = 0


      i = i + 1
      start = 0
      end = 0
   print(distans)
   print(new_row)
   return sum

def matrix_Equals(aoa ,aoa_num, index):
   i = index[0]
   j = index[1]
   new_row = []
   new_row.append([row[0] for row in aoa])
   new_row.append([row[1] for row in aoa])
   new_row.append([row[i] for row in aoa])
   new_row.append([row[j] for row in aoa])
   dis = ArrayEquals([row[i-2] for row in aoa_num], [row[j-2] for row in aoa_num])
   x = expect_roh(dis,new_row)
   return (x)

def selectTWO(aoa,Cow,Sire):
   cow = -2
   sire = -2
   counter = -2
   for row in aoa[0]:
      if row != "":
         counter += 1
      if row == str(Cow):
         cow = counter
      elif row == str(Sire):
         sire = counter
   return cow,sire



def madadim():
   df = pd.read_excel('Moran_cows_if.xlsx')
   new_roh = []
   x = new_matrix(aoa)

   for i in range(len(df["Cow"])):
      y = selectTWO(aoa_all, df["Dam"][i], df["Sire"][i])
      roh_list = matrix_Equals(aoa_for_roh, x, y)
      new_roh.append(roh_list)
   df["expect-roh"] = new_roh

   new_roh = []
   for i in range(len(df["Cow"])):
      y = location_dam(aoa_all, df["Cow"][i])
      roh_list = roh_individual(aoa_for_roh,x, y)
      new_roh.append(roh_list)
   df["ROH-Cow"] = new_roh

   new_roh = []
   for i in range(len(df["Cow"])):
      y = location_dam(aoa_all, df["Sire"][i])
      roh_list = roh_individual(aoa_for_roh,x, y)
      new_roh.append(roh_list)
   df["ROH-Sire"] = new_roh

   new_roh = []
   for i in range(len(df["Cow"])):
      y = location_dam(aoa_all, df["Dam"][i])
      roh_list = roh_individual(aoa_for_roh,x, y)
      new_roh.append(roh_list)
   df["ROH-DAM"] = new_roh
   df.to_excel("Measures.xlsx")

#madadim()

