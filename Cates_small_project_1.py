#Hazelyn Cates
#EN.605.651.81.FA23
#Started 9/7/23
#This program finds the LOCAL alignment of two DNA sequences using a dynamic programming approach
#and prints the alignment score and the two aligned subsequences

import re


def score(dp_table, match_table, seq1, seq2, seq1_len, seq2_len): #this function fills out the DP table
    print("Please enter gap penalty:")
    gap = int(input())
    #print(gap)

    print("Please enter match score:")
    match = int(input())
    #print(match)

    print("Please enter mismatch score:")
    mismatch = int(input())
    #print(mismatch)

    for i in range(1, seq1_len + 1):
        for j in range(1, seq2_len + 1):
            #Checking first condition of Mi-1,j-1 + s(aij), which calculates val1:
            if seq1[i - 1] == seq2[j - 1]:
                val1 = dp_table[i - 1][j - 1] + match

            else:  #meaning it's a mismatch
                val1 = dp_table[i - 1][j - 1] + mismatch

            #print(val1)

            #values 2 and 3 are irrelevant to if the sequences match or not
            #And once all 3 values are calculated, the index i,j in the matrix gets populated w/ the max of the 3
            #Or zero if none of the values are positive (i.e. the max is still a negative number)

            #Checking second condition of Mi-1,j - gap:
            val2 = dp_table[i - 1][j] + gap
            #print(val2)

            val3 = dp_table[i][j - 1] + gap
            #print(val3)

            max_val = max(val1, val2, val3)

            if max_val < 0:  #if the largest number of the 3 is still negative, then that cell in the matrix gets set = 0
                max_val = 0

            print(max_val)

            dp_table[i][j] = max_val

            if max_val == dp_table[i - 1][j - 1] + match: #meaning it is a match
                match_table[i][j] == "diag"
            elif max_val == val2: #Mi-1,j + gap
                match_table[i][j] == "left"
            elif max_val == val3: #Mi,j-1 + gap
                match_table[i][j] == "up"
            else:
                match_table[i][j] == 0


    print(match_table, "\n")

    #Therefore, when doing the traceback, you're only looking at celels that have a "4", meaning the sequences match at that index

    print(dp_table)

    #find max value in DP table
    final_max_val = max(max(i) for i in dp_table)
    print("Max value in DP table = %d" % final_max_val)

    #Extract the index this max value occurs at:
    for i in range(1, seq1_len + 1):
        for j in range(1, seq2_len + 1):
            if dp_table[i][j] == final_max_val:
                x = i #stores index of sequence 1
                y = j #stores  index of sequence 2

    #print("Index of max value: (%d, %d)" % (x, y))

    #So now that I have this max value and it's index, I can write another function that does the traceback and
    #call it in this function
    #Yee haw
    traceback(dp_table, x, y, seq1, seq2)


def traceback(dp_table, x, y, seq1, seq2): #this function traces back the aligned sequence starting from the index of the max value
    #to do this, I need look backwards
    #so like, I need to look at the cells at i-1,j-1, i-1,j and i,j-1 from the cell that has the max value
    #and pick the max of those 3
    #And once one of the 3 values is 0, end at that square and DON'T count that base in the alignment (if it comes before index (0,0))


    """for i in range(x, 0, -1): #don't want to start beyond the index of the max value
        for j in range(y, 0, -1): #same deal, only want to move backwards
            #Also, want to extract the base at that index in the sequence and store in an array or smthn
            #BUT, recall that the sequence indices will be x-1 and y-1
            #So if the max value is at index (2,2) in the DP table,
            #then the matching base will occur at index (1,1) in the sequence
            #yeah, that's right. good god

            #print(i, j)
            #print(max(dp_table[i-1][j-1], dp_table[i-1][j], dp_table[i][j-1]))"""





#First, get the user input:
print("Please enter first sequence:")
seq1 = input()
seq1 = seq1.upper()
#print(seq1)

#Check to make sure that sequence 1 is valid (i.e. only contain A, T, C and G)
if re.search(r'[^ATCG]', seq1):  # checks if any characters BESIDES A, T, C or G are present in sequence 1
    print("Input sequence 1 is not valid. Please restart with valid sequence.")
    exit(1)

print("Please enter second sequence:")
seq2 = input()
seq2 = seq2.upper()
#print(seq2)

#Check to make sure that sequence 2 is valid:
if re.search(r'[^ATCG]', seq2):  # checks if any characters BESIDES A, T, C or G are present in sequence 1
    print("Input sequence 2 is not valid. Please restart with valid sequence.")
    exit(1)

#Find length of each sequence:
seq1_len = len(seq1)
seq2_len = len(seq2)
#print(seq1_len, seq2_len)

#Create DP table that has seq1_len + 1 columns and seq2_len + 1 rows, initialize all cells with zeros
align_table = [[0 for i in range(seq1_len + 1)] for j in range(seq2_len + 1)]  #initialize table to all zeros
#print(align_table)

match_table = [[0 for i in range(seq1_len + 1)] for j in range(seq2_len + 1)]

#Call "score" function to fill out the DP table:
score(align_table, match_table, seq1, seq2, seq1_len, seq2_len)

#Now I need to write a function that does the traceback, but I think I figured that out...
#it seems to be


#Then, I need to somehow implement the arrow thing?? I think that's in the lecture? Lecture 2D maybe?
#That might not be necessary for the traceback, but more for visual purposes when actually drawing the table out
#I dunno, I'll figure it out somehow

#And then I need to print the score and the alignment of the sequences
#Which could perhaps prove to be more difficult than it sounds lol
