import sys
import pandas as pd
import numpy as np
import math
import symnmfModule
np.random.seed(0)
N = -1


def get_goal_type(goal):
    if goal == "symnmf":
        return "SYMNMF"
    elif goal == "sym":
        return "SYM"
    elif goal == "ddg":
        return "DDG"
    elif goal == "norm":
        return "NORM"
    else:
        return "INVALID"


def check_clusters(k, goal, inputFileData):
    global N
    goal_type = get_goal_type(goal)

    if goal_type == "INVALID":
        print("An Error Has Occurred")
        sys.exit(1)

    try:
        with open(inputFileData, "r") as file:
            lines = file.readlines()
            n = len(lines)  
            d = len(lines[0].split(","))
            vectorX = []

            for line in lines:
                values = [float(x) for x in line.strip().split(",")]
                vectorX.append(values)

    except FileNotFoundError:
        print("An Error Has Occurred")
        sys.exit(1)


    N = n
    if (k > 1 and k < N) == False:
        print("An Error Has Occurred")
        sys.exit()

    if goal_type == "SYMNMF":
        W = symnmfModule.norm(vectorX)
        mean=np.mean(W)
        inital_H=np.random.uniform(0,2*math.sqrt(mean/k),size=(n,k))
        inital_H_list = inital_H.tolist()
        return(symnmfModule.symnmf(vectorX,inital_H_list,k, n, d))
 

    if goal_type == "SYM":
        return(symnmfModule.sym(vectorX))
       

    if goal_type == "DDG":
        return(symnmfModule.ddg(vectorX))

    if goal_type == "NORM":
        return(symnmfModule.norm(vectorX))




def main(argv):
    if len(argv) == 4:
        if argv[1].isnumeric() == False:
            print("An Error Has Occurred")
            sys.exit()

      
        if argv[2].isnumeric() == True or (
            argv[2] != "symnmf"
            and argv[2] != "sym"
            and argv[2] != "ddg"
            and argv[2] != "norm"
        ):
            print("An Error Has Occurred")
            sys.exit()
        else:
            k = (int)(argv[1])
            goal = argv[2]
            inputFileData = argv[3]
            answer = check_clusters(k, goal, inputFileData)
            for row in answer:
                print(",".join(str("{:.4f}".format(round(x, 4))) for x in row))


    else:
        print("An Error Has Occurred")
        sys.exit()


if __name__ == "__main__":
    main(sys.argv)





