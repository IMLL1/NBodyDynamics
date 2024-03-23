numbod = 25         # Number of planetary bodies
numdim = 3          # Number of dimensions total
expanded = False    # Whether to expand summations in the dynamics matrix (don't do this for numbod > around 10)
makeTxt = False     # Whether to putput txt file
makeTex = True      # whether to putput Tex file

N = 2*numdim
L = ["x", "y", "z"]

derivOp = "\\frac{\operatorname{d}}{\operatorname{d}t}"
stateVec = "\\begin{bmatrix}"
for b in range(numbod):
    for d in range(numdim):
        stateVec += str(L[d])+"_{"+str(b+1)+"}\\\\"
    for d in range(numdim):
        stateVec += "\\dot{"+str(L[d])+"_{"+str(b+1)+"}}\\\\"
stateVec += "\\end{bmatrix}"

## DYNAMICS MATRIX
dynMtx = "\\begin{bmatrix}\n"
cells = [["0"]*N*numbod for a in range(N*numbod)]
# enter 1s
for bod in range(numbod):
    for dim in range(numdim):
        cells[bod*N+dim][bod*N+dim+numdim] = "1"
# enter off-body gravity terms
for targetnum in range(numbod):
    for sourcenum in range(targetnum):
        for dim in range(numdim):
            cells[targetnum*N+numdim+dim][sourcenum*N+dim] = "\\frac{Gm_{"+str(sourcenum+1)+"}}{r_{"+str(targetnum+1)+","+str(sourcenum+1)+"}^3}"
            cells[sourcenum*N+numdim+dim][targetnum*N+dim] = "\\frac{Gm_{"+str(targetnum+1)+"}}{r_{"+str(sourcenum+1)+","+str(targetnum+1)+"}^3}"
# add same-body gravity terms
for bod in range(numbod):
    for dim in range(numdim):
        if expanded:
            cells[bod*N+dim+numdim][bod*N+dim] = ""
            for sourcenum in range(numbod):
                if bod != sourcenum: cells[bod*N+dim+numdim][bod*N+dim] += "-\\frac{Gm_{"+str(sourcenum+1)+"}}{r_{"+str(bod+1)+","+str(sourcenum+1)+"}^3}"
        else:
            cells[bod*N+dim+numdim][bod*N+dim] = "\\displaystyle\\sum_{n=1,n\\ne"+str(bod+1)+"}^{"+str(numbod)+"}\\frac{-Gm_n}{r_{"+str(bod+1)+",n}^3}"
            
# construct the mtx
for row in range(N*numbod):
    for col in range(N*numbod):
        dynMtx += cells[row][col]
        if(col < N*numbod-1):
            dynMtx += " & "
        else:
            dynMtx += " \\\\ \n"
dynMtx += "\\end{bmatrix}"

# Construct first-order DEs
DEs = "\\begin{aligned}\n"
for b in range(numbod):
    for d in range(numdim):
        DEs += derivOp+str(L[d])+"_{"+str(b+1)+"}&=\\dot{"+str(L[d])+"_{"+str(b+1)+"}} \\\\\n"
        DEs += derivOp+"\dot{"+str(L[d])+"_{"+str(b+1)+"}}&="
        for sourcenum in range(numbod):
            if b != sourcenum: DEs += "-\\frac{Gm_{"+str(sourcenum+1)+"}}{r_{"+str(b+1)+","+str(sourcenum+1)+"}^3}(x_{"+str(b+1)+"}-x_{"+str(sourcenum+1)+"})"
        DEs += "\\\\ \\\\ \n"
DEs += "\\end{aligned}"

# construct second-orer vector DEs
vecDEs = "\\begin{aligned}\n"
for b in range(numbod):
    vecDEs += "\\mathbf{r}_{"+str(b+1)+"}''&="
    for sourcenum in range(numbod):
        if b != sourcenum: vecDEs += "-\\frac{Gm_{"+str(sourcenum+1)+"}\\left(\\mathbf{r}_{"+str(b+1)+"}-\\mathbf{r}_{"+str(sourcenum+1)+"}\\right)}{\\left|\\mathbf{r}_{"+str(b+1)+"}-\\mathbf{r}_{"+str(sourcenum+1)+"}\\right|^3}"
    vecDEs += "\\\\ \\\\ \n"
vecDEs += "\\end{aligned}"
    


## OUTPUT
out = derivOp+stateVec+"=\n"+dynMtx+"\n"+stateVec

# output simple txt
if makeTxt:
    f = open("Dynamics.txt", "w")
    f.write("$$"+out+"$$")
    f.close()

# Construct full LaTeX document
if makeTex:
    header = "\\documentclass{article}\n"
    header += "\\usepackage[margin=25pt]{geometry}\n"
    header += "\\usepackage{amsmath}\n"
    if expanded:
        header += "\\geometry{paperwidth="+str(min(10+1.25*numdim*numbod**2,575))+"cm, paperheight="+str(min(30+1.5*numbod*numdim,575))+"cm}\n"
    else:
        header += "\\geometry{paperwidth="+str(min(10+4*numdim*numbod,575))+"cm, paperheight="+str(min(30+1.5*numbod*numdim,575))+"cm}\n"
    header += "\\title{System Dynaics for "+str(numbod)+"-Body System in "+str(numdim)+" Dimensions with Newtonian Point-Masses}\n"
    header += "\\author{Polaris via Python Code}\n"
    header += "\\begin{document}\n"
    header += "\\setlength{\\arraycolsep}{5pt}\n"
    header += "\\setcounter{MaxMatrixCols}{+"+str(N*numbod+1)+"}\n\n"
    header += "\\maketitle\n\n"
    footer = "\n\\end{document}"

    fullDoc = header
    fullDoc += "\\section*{\\centering{Linearized State Space Model}}\n"
    fullDoc += "\\["+out+"\\]"
    fullDoc += "\\newpage\n"
    fullDoc += "\n\\section*{\\centering{First Order ODEs}}\n"
    fullDoc += "\\["+DEs+"\\]"
    fullDoc += "\\newpage\n"
    fullDoc += "\n\\section*{\\centering{Vector Form}}\n"
    fullDoc += "\\["+vecDEs+"\\]"
    fullDoc += footer

    f = open("Dynamics.tex", "w")
    f.write(fullDoc)
    f.close()
