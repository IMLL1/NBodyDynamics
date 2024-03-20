numbod = 10
numdim = 3
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
            cells[targetnum*N+numdim+dim][sourcenum*N+dim] = "\\frac{Gm_{"+str(sourcenum+1)+"}}{r_{"+str(targetnum+1)+","+str(sourcenum+1)+"}^2}"
            cells[sourcenum*N+numdim+dim][targetnum*N+dim] = "\\frac{Gm_{"+str(targetnum+1)+"}}{r_{"+str(sourcenum+1)+","+str(targetnum+1)+"}^2}"
# add same-body gravity terms
for bod in range(numbod):
    for dim in range(numdim):
        cells[bod*N+dim+numdim][bod*N+dim] = ""
        for sourcenum in range(numbod):
            if bod != sourcenum: cells[bod*N+dim+numdim][bod*N+dim] += "-\\frac{Gm_{"+str(sourcenum+1)+"}}{r_{"+str(bod+1)+","+str(sourcenum+1)+"}^2}"
# construct the mtx
for row in range(N*numbod):
    for col in range(N*numbod):
        dynMtx += cells[row][col]
        if(col < N*numbod-1):
            dynMtx += " & "
        else:
            dynMtx += " \\\\ \n"
dynMtx += "\\end{bmatrix}"


## OUTPUT
out = "$$"+derivOp+stateVec+"=\n"+dynMtx+"\n"+stateVec+"$$"
header = "\\documentclass{article}\n"
header += "\\usepackage[margin=0pt]{geometry}\n"
header += "\\usepackage{amsmath}\n"
header += "\\thispagestyle{empty}\n"
header += "\\pdfpageheight"+str(15 + numbod)+"cm\n"
header += "\\pdfpagewidth"+str(15 + 1.25*numbod)+"cm\n"
header += "\\title{System dynaics for "+str(numbod)+"-Body system in "+str(numdim)+" Dimensions}\n"
header += "\\author{Polaris's Code}\n"
header += "\\begin{document}\n\\maketitle\n\n"
footer = "\n\end{document}"

f = open("Dynamics.text", "w")
f.write(out)
f.close()

f = open("Dynamics.tex", "w")
f.write(header+out+footer)
f.close()