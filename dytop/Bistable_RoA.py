# Bistable_RoA.py  # 2022-20-01
# MIT LICENSE 2020 Ewerton R. Vieira
"""Convert any system with two attractors into a bistable MG"""
import RoA
import matplotlib.pyplot as plt

def assignment(order_retraction):
    assign = dict()
    for k, v in order_retraction.items():
        for i in v:
            assign[i] = k
    assign[-1]=-1
    return assign

def update_file(order_retraction, file_input, name="bistable_RoA.csv"):
    assign = assignment(order_retraction)

    with open(name, "w") as file_out:
        with open(file_input, "r") as file:
            # for i in range(3):
            #     line = file.readline()
            #     file_out.writelines(line)
            line = file.readline()
            while line != "":
                line = line.split(",")
                if line[0][0:4]!="Tile" and len(line)>3:
                    line[1] = str(assign[int(line[1])])
                line = ",".join(line)
                file_out.writelines(line)
                line = file.readline()



order_retraction = {1:{3}, 0:{0,1,2,8,5,6}, 2:{4,7,9,10,11}}
name = "output/bistable_RoA_.csv"

assign = assignment(order_retraction)

update_file(order_retraction, "output/MG_RoA_.csv", name)



lower_bounds = [-1,-1]
upper_bounds = [1,1]
fig, ax = RoA.PlotTiles(lower_bounds, upper_bounds,
            from_file="bistable")

out_pic = "output/bistable_RoA_"

plt.savefig(out_pic, bbox_inches='tight')