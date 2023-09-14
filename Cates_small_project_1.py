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

    for i in range(1, seq2_len + 1): #rows
        for j in range(1, seq1_len + 1): #columns
            #print(seq1[j-1], seq2[i-1])
            #Checking first condition of Mi-1,j-1 + s(aij), which calculates val1:
            if seq2[i - 1] == seq1[j - 1]:
                val1 = dp_table[i - 1][j - 1] + match
                #print("There is a match at index %d in seq2 and index %d in seq1" % (i, j))

            else:  #meaning it's a mismatch
                val1 = dp_table[i - 1][j - 1] + mismatch

            #print(val1)

            #values 2 and 3 are irrelevant to if the sequences match or not
            #And once all 3 values are calculated, the index i,j in the matrix gets populated w/ the max of the 3
            #Or zero if none of the values are positive (i.e. the max is still a negative number)

            #FIX THIS!!!! You know what, I think the issue here is my indices
            #Like, my i and j are switched
            #B/c for example, row 2, column 3 has a max value of six and should be a left arrow
            #where the value comes from the same row but the previous column
            #BUT, with the way I have my loop set up, j = seq1 (columns) and i = seq2 (rows)
            #So it's being interpreted as the opposite, so i should be j and j should be i

            #Checking second condition of Mi-1,j - gap:
            val2 = dp_table[i][j - 1] + gap

            val3 = dp_table[i - 1][j] + gap

            max_val = max(val1, val2, val3)

            if max_val < 0:  #if the largest number of the 3 is still negative, then that cell in the matrix gets set = 0
                max_val = 0

            #print("Max value of row %d of seq2, column %d of seq1 = %d" % (i, j, max_val))

            dp_table[i][j] = max_val #all the numbers are correct, it's "match_table" that's wrong...

            if seq2[i - 1] == seq1[j - 1]: #meaning it is a match
                match_table[i][j] = "diag"
            elif max_val == val2: #Mi-1,j + gap
                match_table[i][j] = "left"
            elif max_val == val3: #Mi,j-1 + gap
                match_table[i][j] = "up"
            else: #mismatch
                match_table[i][j] = 0

    print(match_table, "\n")

    #Therefore, when doing the traceback, you're only looking at cells that have "diag", meaning the sequences match at that index

    print(dp_table)

    #find max value in DP table
    final_max_val = max(max(i) for i in dp_table)
    print("Max value in DP table = %d" % final_max_val)

    #Extract the index this max value occurs at:
    for i in range(1, seq2_len + 1):
        for j in range(1, seq1_len + 1):
            if dp_table[i][j] == final_max_val:
                x = i #stores index of sequence 2
                y = j #stores  index of sequence 1

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

    #Character list to hold to aligned bases:
    aligned_seq = list()
    #Lists to keep track of the indices where the matches occurred in each sequence
    seq1_bases = list()
    seq2_bases = list()

    #keep track of the alignment score:
    score = 0

    while dp_table[x][y] != 0:
        if match_table[x][y] == "diag":
            seq1_bases.append(y - 1)
            seq2_bases.append(x - 1)

            score += match_table[x][y]
            x -= 1
            y -= 1

    """for i in range(x, 0, -1): #don't want to start beyond the index of the max value
        for j in range(y, 0, -1): #same deal, only want to move backwards
            if match_table[i][j] == "diag":
                #Record the index the match occurred in each sequence:
                seq1_bases.append(j-1)
                seq2_bases.append(i-1)
                #Make sure to append any aligning base to the FRONT of the list since the traceback starts at the end

                #FIX THIS!!!! It's not appending the correct bases!!
                aligned_seq[:0] = seq1[i - 1]  #since the bases match, it doesn't matter which sequence you reference for the nucleotide
                i = i - 1 #move up one row
                j = j - 1 #move left one column and start from there
                print(i, j) #FIX THIS!!! Why is it not decrementing j???

                if dp_table == 0: #meaning it's the end of the alignment
                    break
                else:
                    score += dp_table[i][j]

                #Now in this new index, check the arrow
                #NOTE: only

            elif match_table[i][j] == "up": #NOT a match, stay in same column, go to previous row
                #decrement j and keep going
                if dp_table == 0: #meaning it's the end of the alignment
                    break
                else:
                    score += dp_table[i][j]
                    i -= 1
                    continue

            elif match_table[i][j] == "left": #NOT a match, go to previous column, stay in same row
                if dp_table == 0: #meaning it's the end of the alignment
                    break
                else:
                    score += dp_table[i][j]
                    #decrement i
                    j -= 1
                    continue

            elif match_table[i][j] == 0:
                if dp_table == 0: #meaning it's the end of the alignment
                    break
                else:
                    i = i - 1
                    j = j - 1
                    score +=  dp_table[i][j]
                    continue"""

                #IMPORTANT!!!!!
                #PROBLEM: for each run through the loop, there can only be ONE match.
                #Meaning, if sequence 2 (vertical) has an A as the last base (since it works from the bottom up)
                #Then it can only match to ONE 'A' in sequence 1 (horizontal) b/c you CAN'T match one base to two bases in an alignment
                #So a "diag" at j = 5 can only match the FIRST (technically the last) "A" it finds in sequence 1.
                #So once a match is found, it has to decrement to the next column even if it didn't get through all the rows
                #UPDATE: RESOLVED: added a "break" to the inner loop


            #Also, want to extract the base at that index in the sequence and store in an array or smthn
            #BUT, recall that the sequence indices will be x-1 and y-1
            #So if the max value is at index (2,2) in the DP table,
            #then the matching base will occur at index (1,1) in the sequence
            #yeah, that's right. good god

            #print(i, j)
            #print(max(dp_table[i-1][j-1], dp_table[i-1][j], dp_table[i][j-1]))

    #These tell you which indices of sequences 2 and 1 (respectively) the matches occur at
    print(seq1_bases)
    print(seq2_bases)

    print("Alignment score = %d" % score)

    print("\nOriginal sequences:")
    print("S1: ", seq1)
    print("S2: ", seq2)

    print("The aligned sequence is:")
    print(aligned_seq)

    """print("Occurring at indices ", end="")
    print(sorted(seq1_bases), " in sequence 1 and indices ", end="")
    print(sorted(seq2_bases), " in sequence 2.")"""


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
    exit(1) #if sequence is invalid, exit program and start again

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
