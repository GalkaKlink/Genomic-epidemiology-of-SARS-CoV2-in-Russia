import sys
import ete3
import PyQt5
from ete3 import Tree, NodeStyle, AttrFace, faces, TreeStyle, SeqMotifFace
import re


with open("Nodes.txt") as f:
    roots = [line.split()[0] for line in f] 
f.close()

def layout(node):
	if node.is_leaf():
		
		if ("Russia/" in node.name) or ("FMBA" in node.name):
			Rudesc = faces.AttrFace("name", fsize=10)   #Dynamic text Face. Text rendered is taken from the value of a given node attribute
			Rudesc.margin_left = 5
			Rudesc.margin_right = 5
			Rudesc.background.color = "#e77576"
			
			faces.add_face_to_node(Rudesc, node, 0, aligned=False)
		
		else:
			if "_" in node.name:
				node.name = node.name.replace("_"," ")
			countrydesc = faces.AttrFace("name", fsize=10)
			countrydesc.margin_left = 5
			countrydesc.margin_bottom = 10
			countrydesc.margin_top = 5
			countrydesc.margin_right = 5
			faces.add_face_to_node(countrydesc, node, 0, aligned = False)


#color = "#3690c0"#input("please enter color:")
color = "black"
#########################
##visualize
ts = TreeStyle()
ts.branch_vertical_margin = 0.5
ts.root_opening_factor = 1
ts.scale =  50
ts.draw_guiding_lines =True   #Draw guidelines from leaf nodes to aligned faces
ts.guiding_lines_type = 2 #dotted


ts.show_leaf_name = False
ts.layout_fn = layout
nstyle = NodeStyle()
#rnstyle = NodeStyle()
lnstyle = NodeStyle()
#dlstyle = NodeStyle()
#color = random_color()
for i in roots:
    tree1 = ete3.Tree(i+".nwk",format=1)
    for rn in tree1.traverse():
        nstyle["vt_line_color"] = color
        nstyle["hz_line_color"] = color
        nstyle["vt_line_width"] = 4
        nstyle["hz_line_width"] = 4
        nstyle["size"] = 0
        rn.set_style(nstyle)
        if rn.is_leaf():
            lnstyle["vt_line_color"] = color
            lnstyle["hz_line_color"] = color
            lnstyle["vt_line_width"] = 4
            lnstyle["hz_line_width"] = 4
            lnstyle["size"] = 5
            lnstyle["fgcolor"] = "black"
            rn.set_style(lnstyle)

    
    outfile = i+".red.pdf"
    tree1.render(outfile, w = 4000, units= 'px', dpi = 300, tree_style = ts)
