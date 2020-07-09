with open("side.txt", "r") as sidestream:
    with open("layer.txt", "r") as layerstream:
        with open("column.txt", "r") as columnstream:
            with open("answers.txt", "w") as outstream:
                for sideline, layerline, columnline in zip(sidestream, layerstream, columnstream):
                    sidelist = sideline.split(",")
                    layerlist = layerline.split(",")
                    columnlist = columnline.split(",")
                    for s, l, c in zip (sidelist, layerlist, columnlist):
                        coord = str(s) + "    " + str(l) + "    " + str(c)
                        print(coord)
                        outstream.write(coord)
                        outstream.write("\n")
