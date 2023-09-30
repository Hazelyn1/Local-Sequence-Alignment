#Hazelyn Cates
#EN.605.651.81.FA23
#Started 9/7/23, finished 9/20/23
#This program finds the local alignment of two DNA sequences using a dynamic programming approach
#and prints the alignment score and the two aligned subsequences

import re

#this function fills out the DP table
def score(dp_table, match_table, seq1, seq2, seq1_len, seq2_len):
    print("Please enter gap penalty:")
    gap = int(input())
    #print(gap)

    print("Please enter match score:")
    match = int(input())
    #print(match)

    print("Please enter mismatch score:")
    mismatch = int(input())
    #print(mismatch)

    for i in range(1, seq2_len + 1): #rows
        for j in range(1, seq1_len + 1): #columns
            #Checking first condition of Mi-1,j-1 + s(aij), which calculates val1:
            if seq2[i - 1] == seq1[j - 1]: #meaning it's a match
                val1 = dp_table[i - 1][j - 1] + match

            else:  #meaning it's a mismatch
                val1 = dp_table[i - 1][j - 1] + mismatch

            #Checking second condition of Mi-1,j - gap:
            val2 = dp_table[i][j - 1] + gap

            val3 = dp_table[i - 1][j] + gap

            max_val = max(val1, val2, val3)

            if max_val < 0:  #if the largest number of the 3 is still negative, then that cell in the matrix gets set = 0
                max_val = 0

            dp_table[i][j] = max_val #set that index of the DP table to the calculated max value

            #Populate the match table (which contains the arrows) according to the max value above
            if seq2[i - 1] == seq1[j - 1]: #meaning it is a match
                match_table[i][j] = "diag"
            elif max_val == val2: #Mi-1,j + gap
                match_table[i][j] = "left"
            elif max_val == val3: #Mi,j-1 + gap
                match_table[i][j] = "up"
            else: #mismatch
                match_table[i][j] = 0

    #find max value in DP table
    final_max_val = max(max(i) for i in dp_table)

    #Extract the index this max value occurs at:
    for i in range(1, seq2_len + 1):
        for j in range(1, seq1_len + 1):
            if dp_table[i][j] == final_max_val:
                x = i #stores index of sequence 2 (rows)
                y = j #stores  index of sequence 1 (columns)

    #Call the traceback function, which finds the aligned subsequence
    traceback(dp_table, x, y, seq1, seq2, final_max_val)


#this function traces back the aligned sequence starting from the index of the max value in the DP table
#And uses the arrows in the match table to generate the local alignment
def traceback(dp_table, x, y, seq1, seq2, max_val):
    #Character list to hold to aligned bases of each sequence:
    aligned_seq1 = list()
    aligned_seq2 = list()

    #keep track of the alignment score:
    score = 0

    while dp_table[x][y] != 0: #iterate through DP table until a score of 0 is hit
        if match_table[x][y] == "diag": #meaning it's a match
            aligned_seq1[:0] = seq1[y-1] #add the base to subsequence 1
            aligned_seq2[:0] = seq2[x-1] #add the base to subsequence 2

            #score += dp_table[x][y] #update the score according to the current value in the DP table at the given index
            x -= 1 #move up one row
            y -= 1 #move left one column

        elif match_table[x][y] == "up": #same column, one row up
            #score += dp_table[x][y] #update score
            aligned_seq1[:0] = "_" #add a gap to subsequence 1
            aligned_seq2[:0] = seq2[x-1] #add base to subsequence 2
            x -= 1 #go up one row

        elif match_table[x][y] == "left": #same row, one column left
            #score += dp_table[x][y] #update score
            aligned_seq1[:0] = seq1[y-1] #add base to subsequence 1
            aligned_seq2[:0] = "_" #add gap to subsequence 1
            y -= 1 #go left one column

        else: #when there's a mismatch and match_table[x][y] == 0
            #score += dp_table[x][y] #update score
            aligned_seq1[:0] = seq1[y-1] #add base to subsequence 1
            aligned_seq2[:0] = seq2[x-1] #add base to subsequence 2
            x -= 1 #go up one row
            y -= 1 #go left one column

    #Print results:
    print("\nAlignment score = %d" % final_max_val)

    print("\nOriginal sequences:")
    print("S1: ", seq1)
    print("S2: ", seq2)

    print("\nAligned subsequences:")
    print("Subsequence 1: ", end="")
    print(aligned_seq1)
    print("Subsequence 2: ", end="")
    print(aligned_seq2)



#First, get the user input:
print("Please enter first sequence:")
seq1 = input()
seq1 = seq1.upper()

#Check to make sure that sequence 1 is valid (i.e. only contain A, T, C and G)
if re.search(r'[^ATCG]', seq1):  #checks if any characters BESIDES A, T, C or G are present in sequence 1
    print("Input sequence 1 is not valid. Please restart with valid sequence.")
    exit(1)

print("Please enter second sequence:")
seq2 = input()
seq2 = seq2.upper()

#Check to make sure that sequence 2 is valid:
if re.search(r'[^ATCG]', seq2):  # checks if any characters BESIDES A, T, C or G are present in sequence 1
    print("Input sequence 2 is not valid. Please restart with valid sequence.")
    exit(1) #if sequence is invalid, exit program and start again

#Find length of each sequence:
seq1_len = len(seq1)
seq2_len = len(seq2)

#Create DP table that has seq1_len + 1 columns and seq2_len + 1 rows, initialize all cells with zeros
align_table = [[0 for i in range(seq1_len + 1)] for j in range(seq2_len + 1)]  #initialize table to all zeros

#Create match table, which will hold to the arrow assignments, which is the same size as the DP table
match_table = [[0 for i in range(seq1_len + 1)] for j in range(seq2_len + 1)] #initialize to all zeros

#Call "score" function to fill out the DP table:
score(align_table, match_table, seq1, seq2, seq1_len, seq2_len)
