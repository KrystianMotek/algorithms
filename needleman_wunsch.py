import numpy as np


class SequenceAlignment:
    def __init__(self, x, y, match, mismatch, indel):
        self.x = list(x)
        self.y = list(y)
        self.match = match 
        self.mismatch = mismatch
        self.indel = indel 

    def build_table(self):
        rows = len(self.y) + 1
        columns = len(self.x) + 1
        table = np.zeros(shape=(rows, columns))

        for i in range(rows):
            table[i][0] = -i 
        
        for j in range(columns):
            table[0][j] = -j

        for i in range(len(self.y)):
            for j in range(len(self.x)):
                row_index = i + 1
                column_index = j + 1
                if self.y[i] == self.x[j]:
                    # match
                    top = table[row_index-1][column_index] + self.indel
                    down = table[row_index][column_index-1] + self.indel
                    diagonal = table[row_index-1][column_index-1] + self.match
                else:
                    # mismatch
                    top = table[row_index-1][column_index] + self.indel
                    down = table[row_index][column_index-1] + self.indel
                    diagonal = table[row_index-1][column_index-1] + self.mismatch
                    
                table[row_index][column_index] = np.max([top, down, diagonal])
        
        return table
    
    def align_sequences(self):
        table = self.build_table()
        i = len(self.y)
        j = len(self.x)
        x_result = ""
        y_result = ""
        while i > 0 or j > 0:
            if i == 1 and j == 0:
                x_result += "-"
                y_result += self.y[i-1]
                i -= 1
            if i == 0 and j == 1:
                x_result += self.x[j-1]
                y_result += "-"
                j -= 1
            else:
                if table[i-1][j-1] >= table[i-1][j] and table[i-1][j-1] >= table[i][j-1]:
                    x_result += self.x[j-1]
                    y_result += self.y[i-1]
                    i -= 1
                    j -= 1
                    continue
                if table[i-1][j] > table[i-1][j-1] and table[i-1][j] >= table[i][j-1]:
                    x_result += "-"
                    y_result += self.y[i-1]
                    i -= 1
                    continue
                if table[i][j-1] > table[i-1][j-1] and table[i][j-1] > table[i-1][j]:
                    x_result += self.x[j-1]
                    y_result += "-"
                    j -= 1
                    continue
                    
        x_aligned = "".join(reversed(x_result))
        y_aligned = "".join(reversed(y_result))
 
        return x_aligned, y_aligned
    

if __name__ == "__main__":
    x = "GCATGCG"
    y = "GATTACA"
    alignment = SequenceAlignment(x, y, match=1, mismatch=-1, indel=-1)
    output = alignment.align_sequences()
    x_aligned = output[0]
    y_aligned = output[1]
    print("Aligned sequences")
    print(x_aligned)
    print(y_aligned)
