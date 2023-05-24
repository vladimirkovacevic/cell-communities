import base64
import os

from .utils import timeit

def get_base64_encoded_img(path):
    data_uri = base64.b64encode(open(path, 'rb').read()).decode('utf-8')
    return '<img src="data:image/png;base64,{0}">'.format(data_uri)

def get_css():
    return r'''
        img { width:100%; height:100%}
        table, th, td {
        border: 0px solid black;
        }
        .twobytwo {
            border: 1px solid black;
        }
        hr {
        width: 100%;
        }
        .rightpad {
            padding: 3px 15px;
        }
        .shrink > img {
            max-width: 50vw;
        }
        .shrinkmore > img {
            max-width: 40vw;
        }
        .tableSmaller img {
            max-width: 70vw;
        }
        .lastTable img {
            max-width: 60vw;
        }
        table.center {
            margin-left: auto; 
            margin-right: auto;
        }
    '''

def get_slices_content(path, plotting_level):
    content = ""
    for root, dirs, files in os.walk(path): 
        for name in dirs:
            content += get_slice_plots(os.path.join(root, name))
            content += get_table_plots(os.path.join(root, name))
            if plotting_level > 2:
                content += get_cmixtures(os.path.join(root, name))
                content += get_colorplots_rgb(os.path.join(root, name))
    return content

def get_slice_plots(path):
    for root, dirs, files in os.walk(path):
        for name in files:
            if name.startswith("clusters_cellspots"):
                clusters = os.path.join(root, name)
            if name.startswith("cell_type_annotation"):
                anno = os.path.join(root, name)
            
    return f'''
        <thead>
            <tr>
                <th style="text-align:left"><h3>Name: {os.path.basename(path)}</h3></th>
            </tr>
        </thead>
        <table>
            <tbody>
            <tr>
                <td>{get_base64_encoded_img(anno)}</td>
                <td>{get_base64_encoded_img(clusters)}</td>
            </tr>
            <tr>
                <td>Plot represents cell types according to the annotation used</td>
                <td>Plot represents detected communities</td>
            </tr>
            </tbody>
        </table>
        <hr>
    '''

def get_table_plots(path):
    for root, dirs, files in os.walk(path):
        for name in files:
            if name.startswith("celltype_table"):
                anno_table = os.path.join(root, name)
            if name.startswith("cell_mixture_table"):
                cell_mixtures = os.path.join(root, name)
    return f'''
        <table>
            <tbody>
            <tr>
                <td class="tableSmaller">{get_base64_encoded_img(anno_table)}</td>
                <td>{get_base64_encoded_img(cell_mixtures)}</td>
            </tr>
            <tr>
                <td>Table represents ...</td>
                <td>Table represents percentages of cell types in each community</td>
            </tr>
            </tbody>
        </table>
        <hr>
    '''

def get_cmixtures(path):
    cmixtures = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if name.startswith("cmixtures"):
                cmixtures.append(os.path.join(root, name))
    tds = []
    for cmix in cmixtures:
        tds.append(f'<td class="shrink twobytwo">{get_base64_encoded_img(cmix)}</td>')

    rows = ""
    for i in range(0,len(tds),2):
        row = "".join(tds[i:i+2])
        rows += f"<tr>{row}</tr>"
    
    return f'''
        <table>
            <tbody>
            <tr>
                <td>These plots represent cell type mixture (left) of obtained communities (right)</td>
            </tr>
            {rows}
            </tbody>
        </table>
        <hr>
    '''

def get_colorplots_rgb(path):
    colorplots = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if name.startswith("colorplot_rgb"):
                colorplots.append(os.path.join(root, name))
    tds = []
    for cp in colorplots:
        tds.append(f'<td class="shrinkmore twobytwo">{get_base64_encoded_img(cp)}</td>')

    rows = ""
    for i in range(0,len(tds),2):
        row = "".join(tds[i:i+2])
        rows += f"<tr>{row}</tr>"

    return f'''
    <table>
        <tbody>
        <tr>
            <td>RGB Plots</td>
        </tr>
        {rows}
        </tbody>
    </table>
    <hr>
    '''

def get_overall_content(path, plotting_level):
    if plotting_level < 3:
        return
    rows = ""
    for item in os.listdir(path):
        if os.path.isfile(os.path.join(path, item)):
            if item.startswith("cell_abundance_all_slices"):
                rows += f"<tr><td>{get_base64_encoded_img(os.path.join(path, item))}</td></tr>"
                rows += "<tr><td>Percentages of each cell type within all slices aligned next to each other for comparison</td></tr>"
            if item.startswith("cell_abundance_per_slice"):
                rows += f"<tr><td>{get_base64_encoded_img(os.path.join(path, item))}</td></tr>"
                rows += "<tr><td>Percentages of each cell type in each slice</td></tr>"
            if item.startswith("cluster_abundance_all_slices"):
                rows += f"<tr><td>{get_base64_encoded_img(os.path.join(path, item))}</td></tr>"
                rows += "<tr><td>Percentages of each community within all slices aligned next to each other for comparison</td></tr>"
            if item.startswith("cluster_abundance_per_slice"):
                rows += f"<tr><td>{get_base64_encoded_img(os.path.join(path, item))}</td></tr>"
                rows += "<tr><td>Percentages of each community in each slice</td></tr>"
            if item.startswith("total_cell_mixtures_table"):
                rows += f"<tr><td>{get_base64_encoded_img(os.path.join(path, item))}</td></tr>"
                rows += "<tr><td>Table represents percentages of cell types in each community for all slices</td></tr>"
            if item.startswith("cell_perc_in_community_per_slice"):
                rows += f'<tr><td>class="lastTable">{get_base64_encoded_img(os.path.join(path, item))}</td></tr>'
                rows += '<tr><td>Table represents total number of cells in each community per slice</td></tr>'
            rows += '<tr style="padding-bottom:50px"><td><br/></td></tr>'
            
    return f'''
    <table>
        <tbody>
        <tr>
            <td>Plots obtained using all slices</td>
        </tr>
        {rows}
        </tbody>
    </table>
    <hr>
    '''

@timeit
def generate_report(args, algo_list):
    if args.plotting == 0:
        return
    
    htmlstr = f'''
    <!DOCTYPE html>

    <html lang=”en”>

    <head>
    <style>
    {get_css()}
    </style>

    <title>Cell Communities Report</title>
    <link rel="icon" type="image/png" sizes="16x16" href="data:image/png;base64,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0AYAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAP///Rbs9P69mMj//4jA//+Owv//ocz//4O9//92tf//g73//6XP//+Qw///frr//3K0//94t///j8P//4nA//+dy///n8v1/5rBnf+iwnj/pcV+/7DLjv+vxIr/r8KJ/7PCjv+qxIT/yLWh/8C5mf/kmLT/3Zut/6mzfP+gwXf/u8mZ/9+qtv/4mMv/27K1/9qor//9i8v/+JvV8ui56WvLq/AGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAP///mvA3f78a7D//3W1//+o0P//ncr//5HE//99uv//d7b//4i///+Dvf//hL3//4G8//+CvP//f7v//4vB//+y1v//h7///5XBt/+dwHD/qsiG/6vGhv/ToaX/vLOR/6PGff+gwnj/p8eC/6DDeP+twIb/tsuV/57Cdf+jxHz/p8eC/6vEhf/Etp3/5afJ//CSzP/4lNL/36Hl/9ms7fPar+2E16fsCQAAAAAAAAAAAAAAAAAAAAAAAAAA///9Dc7l/7+Pw///brH//4vB//+jzf//hL3//5bH//97uf//b7L//4C7//+QxP//fLn//3W1//97uf//Z67//5XG//+jzv//i8H6/7HRxP+0z5H/rcqI/7LMkf/cqLL/06am/7PNkv+wzI7/qsiG/6HCef+uy4v/r8yN/7HHjf+mxH//nsF2/6XGf//AxK//3b/r/+eo5f/zmtf/2ajs/9ux7//euPDn3rjwKgAAAAAAAAAAAAAAAAAAAAAAAAAA3+3+QYzC//FkrP//c7T//4S9//+Ow///g73//4S9//91tf//drb//2qw//9gqv//ebf//3+7//96uP//fbn//4W+//+Hv///kL/I/7DMjf+uzJ3/psie/67Kjv/AsJT/wq6W/6zJif+lxX//sc2Q/67LjP+mxYD/zrqo/+SywP+7ypr/sMWM/9a5s//mqd//3LHv/+Kw6v/opuL/3bLu/9qj6f/qsOb49afbUwAAAAAAAAAAAAAAAAAAAAD///oF3Oz+oqrS//+ey///lMb//4e///+Dvf//h7///4S9//+Iv///iL///2qw//9nrv//b7L//4G8//90tf//ksX//5/M//+RxP//lsbt/53H0P+Buuv/l8TL/6jHhv+pyIX/rcyL/6bGgP+uy4v/rsuL/6fHgv+7wJb/9JzI//Cixv/Kzq3/36u2//yJyf/1kNL/26Hn/9ux7//Xqe3/4KTm//OGzv/7hcn//JHNu/qa0BMAAAAAAAAAAAAAAAD0+P0tqtL/5326//+Gvv//ebf//3y5//+Bu///jsP//3a2//+RxP//jcL//4O9//9ws///e7j//4a+//90tP//kMP//5HE//+Sxf//eLf//3q4//+Iv/f/r82m/6TEfP+1z5b/wc2j/6XBff+YvW3/nsB0/6XFf/+xyY7/sbGF/7y9lf/BvZv/74y8//54xP/6hsf/5pnc/9uq7P/opeL/95XS//eT0v/yrOD//JTP+ft6wlcAAAAAAAAAAAAAAADP5f5mjsL//o/D//+RxP//lcb//5fH//90tf//fLn//4zC//+BvP//ksX//4nA//+Iv///hr7//366//+Wx///jML//4O9//+u1P//jMDr/3u25f+MwPP/tdCg/6bGgP+8zp3/1721/7G/iv+ixXv/nsB0/5/Bdv+lxX//pseB/6HEev+uyov/wLSW/92arf/Rqaf/3LDS/96q6//2ndf//pzS//Gw4//iue//8pfW//mX0qj9wt8IAAAAAP///QW93P+ogrz//5rJ//+cyv//oc3//4vB//+CvP//mMj//6PO//+Au///icD//43C//94t///eLf//4zC//+Dvf//bLD//3y5//+cyv//o8i9/63Nuv+Lwf3/m8fY/7jSov+qx4X/q8SF/6O6eP+7vpX/scmO/6XEfv+bvnD/pcaA/6XFfv+kxH3/nsF1/8LVp/+/06b/xbHA/+ac3v/5kc///YvL//id1v/krOj/3q3r/+Su6OXkuuwqAAAAAP///Ryz1f7VdLT//3+6//+Lwf//grz//3u5//+Hv///mcj//5DE//+Dvf//g73//4vB//+Fvv//erj//3Cy//9ztP//hb7//4rA//9/u///msTH/6zMq/+Bu/f/jLy4/7fRmv+tyov/o8F7/9qpsP/6qNT/7ajG/9S0rv+nxIH/xcOj/7jNmP+nyIL/x8Cj/8rFqf+gxHj/sLaH//aJxP/9gMb//H7F//2PzP/6kM//7KDd/9yt7Pnhs+xQAAAAALva/kOs0v/viL///4W+//+Txf//frr//4K8//+Fvv//ksX//366//+Mwv//grz//4nA//9+uv//grz//366//9ztP//kMP//4nA//+Hv///mMfv/6zMsf+PvsD/l76C/6LDef+ox4P/sMaL/96otP/4ksn//5PR/+Sht/+swob/6qfC/9fBt//BuZr/76HF/9qosP+owYD/rMWG/9S0rf/7kMz//InK//ej2f/0jND/+qPX//ak2v/4ldKPAAAAAIS9/2h6uP/8lsf//3+7//92tv//ebf//366//+Fvv//gLv//3+7//+ey///cbP+/3Wz+v+EvPv/ebj//2+y//9tsf//jcL//4O8//+Au///dbTw/6nKp/+zzpL/oMJ4/6LDev+tyor/t9CZ/7/To//Pt6n//pnS/+imvv/Au5n/8JC+//yk1f/7ptT//aPW//iSyf/pob7/0bGp/9GwqP/6mM7//ajY/+i56//es+7/8rLi//qm2f/3qNux8J7XCL3b/pORxP//jcL//43C//99uv//e7j//3W1//+Hv///crT//4e///+Cvf//gLn6/7K8wv+8w8f/frn9/2iu//99uf//grz//4S9//+Gvv//ib/5/7DQxP+vzJH/l7xq/57AdP+yzZH/vNSg/5/Cdv+qxYT/4qy7//ut2P/4qtP//J3S//2W0P/5s9//9qzd//2V0P/+mNL/+5HM//uHyP/9kM3//ZrS//eo2//jtOz/26Tp/+Kw6v/qtujk+prRQLPV/L2Owv//ksX//57L//+ZyP//isD//3a2//+Kwf//Z67//3W1//+Evf//iL///6jI6v+svc3/e7b6/2uw//9ws///a7D//3a2//+Kwf//ksX+/4S86f+bwI//sMeK/67Ki/+syYn/sMyO/7LDjf/Fs5v/4pqz//2Pzf/9kM7//Y3M//2j1v/2yOv/9b/n//2t2v/9ntT//Y3M//2W0P/8gMb//JDN//2FyP/vmNj/4Lvx/+G57//zodr+/IPGgaLN/86gzP//oMz//4vB//+Dvf//nMr//5jI//9ytP//b7L//4O9//+QxP//nsv//4O9//96uP3/Z67//3Gy/P+Ku/P/d7b//4jA//+Cu/X/pcrC/6DGof+/v57/2J+p/7XKlP+1z5X/s9CU/8e5of/7lc7//YHH//6Y0v/+mNL//ZPP//2q2f/80+v//Mzp//2q2f/9k87//IbJ//2Iyv/9j83/7KXD/+qnwf/5st3/5c7j/++v4f/8odX9/JbPhoi//9CUxv//l8Hy/3uz8/9mrv//grz//6rR//+Qw///brL//3e2//97uP//icD//4C7//+Evf//Za3//43A/P+jyPH/hb7//4S9//+v0ur/utKf/7bQlP+7y5v/sbeG/7XPlv+ryoj/r86P/9S+sv/5icb/947G/+uavv/iqLn/+I3G//2X0f/9qNj//rzh//2j1v/9ltD//ITI//2DyP/8p9b/1sK3/+qiv//8otX/9q7T//ye0//9ntP1/LbeVXy5/96Ow///s8nf/5S75f9ws///ksX//5nI//+Txf//eLf//326//+KwP//grz//2eu//9ztP//bLH//4e///+Sxf//nsv//5jI//+kze3/vtam/7TOlP+80J7/xs2o/8vAqf/Px7D/v8ie/7vLmv/JxKj/ur2T/7O9jP+1z5X/z7Cm//iVy//5n9D/98vg//qg0f/9k8///ZPO//2d0//8kM3/+JbK//yNzP/9nNL//ZbQ//2X0P/9pNbr/LjeOpHE/86Pw///sNX//6PO//9ytP//hb7//43C//+Evf//eLf//4nA//+ey///j8P//326//+Dvf//hL3//43C//+Zyf//grz//4a+9/+my8P/vtWh/8eynv/yq83/9aTN/+28zv/7nND/6py9/8DEnv+205n/uNKa/7DIjf/ErZj/2Zup//mez//UtK3/x8up/9O9sf/zncf//qDV//2s2v/9icv//ZbR//6+4v/9tt7//pLP//yU0P/8o9b4/KPVZYi//ryLwf//k8X//4vB//+Qw///lsf//4e///9xs///gLv//47C//95t///j8P//5DE//+QxP//ncv//4i///+Uxv//t9j//5HD8P+qypj/s82S/9Gyqf/eq7X/2be1/7jHlf/qmbv/7J/A/77PoP+nyIP/tNCV/6vDhf/cn67/9JjH/+ilvv+1wo//x92w/77Xo//Sx7P/85jG//en0f/3oc7/9qPO//uj0//7otP/356y/+aw2P/vrOL2/I7MXpLE/r2JwP//gLv//3u4//+izf//icD//3Gz//+CvP//hL3//4e///96uP//f7r//3i3//9ztP//iL///325//95t///ncv//4q/6f+vzJf/ss2R/6XGf/+ty4v/rMyK/7jKl//xpMj/8rPQ/8HHoP+7wpf/u8WY/7nSm/+wyY3/wMyg/8nPrf+1u6H/2dDP/863qf+wtoX/56S9/8KymP+7w5j/wr6d/8monP/LsqL/psWA/8fGuP/usuTh+ZbPK3a2/md4t//4hL3//325//+Hv///grz//4K8//+Xx///tNf//53L//91tf//g7z//4/D//9+uv//drb//3i3//+Xx///hb7+/4a73v+20aH/rMqJ/6TGfv+wzI7/t8aU/8PGof/LtaP/2Kyv//Ghx//wnsT/2bq2/9i4tP/Yt7P/yMWn/9qyuf/dt+T/7b/u//eMyP/0lsX//JnQ/9ujr/+5uJD/sMKK/7DCiv+1ypf/zruq/+emwf/8qdjT+abVI1un/it3tv/djsP//3u4//9qr///dLX//3S1//+byv//pM7//4rB//+Evf//j8P//5rJ//+Gvv//gLv//43C//+Hv/7/gbjM/7HOn/+505v/yL6j/8uxov/MsaP/7KbE/+60zP/GxKT/vMyc/++dw//9hsr/+ZfN//2Szv/8k87/8pjE//uWzv/1mdb/96Xa//2b0v/9hsn//ZLO//6Tz//3os7/86nL/+2fzf/bveP/6rjl//mY1Pn8g8dfAAAAAIO8/Q6Lwf/Bhr7//3y5//+BvP//ZKz//3i3//+Txf//hr7//4S9//+Au///k8X//5LF//+Au///lsf//43C//+Lvtv/q8qY/6/MjP+tvYX/6ZW5//2Szv/7mdD//qLX//yX0P/WtLH/z9e3//OZx//8l9D//I/M//2Hyf/9l9D//qHW//2Pzf/3kcz//K3a//yi1v/0o9v/95PS//yi1f/9hMj//JHP/+ut5f/etvD/4L3x/+m46rv5c8EWAAAAAEmc+wN5t/6dps///5HE//+Lwf//gLv//3O0//96uP//h7///4rB//97uP//isD//3e2//+Xx///sNX//4zB+v+mybD/ss2O/67Li/+sv4T/zqqj/+uuxv/VrKv/1r20//mf0P/KvKb/xr2i//ySzf/prcL/3aiz//2V0P/9p9f//YjK//iTzf/Hqcz/7qXT//qIzP/ltev/6rTo//uZ0//9k87/9qba/+G78P/asO7/5cfz/+G+8IUAAAAAAAAAAAAAAACRxP5aos3//I/D//+byf//mMj//3+7//+ey///fbr//4G7//+Owv//crT//3Cz//+Iv///eLf//3248f+oyJn/utKc/6/Ljf+zz5P/utSd/67Fif+jw3v/tsOS//Cty/+/vpr/u8iZ//Ouz//jqrv/4rC9//6Z0//9ldD//YvL//yTzv/qqNT/95/S//yKy//2m9f/8L3a//ii0P/+ntT/8Jza/+K78P/muez/6cry/+LC8YUAAAAAAAAAAAAAAACBu/wvj8P+6q3T//+dyv//icD//325//+Evf//hL3//366//+Dvf//bbH//4vB//+Txf//frr//5rI+f+71K//rcqJ/6LBef/Nr6T/yraj/6LEe/+pyIT/utGd/8POpf+1zpb/v8ac/++cw//8mdD//JjP//2Ky//9jcz//ZLO//2h1f/+otb//aXX//2e0//tp8X/vsGZ/8/FsP/xrcz/5rLe/+e66v/5o9j/67Hm8uDA8T8AAAAAAAAAAAAAAAB7svALjsH8rKHN//94t///c7T//26y//9+uv//frr//2iu//9ytP//gLv//3q4//99uf//jMH//5fH/P/N4uD/ttCV/7e0jP/ymcX/yq6g/6zLif+xzpH/q8mH/67Ki/+vv4f/5Km7//yTz//1mNX/87fk//qZ0v/9kM3//YzL//2Tzv/9lM///YLH//iRyP/NuKf/w8Wh/9XLuf/Px7D/1LWy//it2P/9mtL/+qDVv+rK7hEAAAAAAAAAAAAAAAAAAAAApM7+VYS9//d2tv//i8H//3q4//96uP//gbz//4G8//+Iv///l8f//4a+//+Bu///lcb//4/D//+41+T/usiX/+mXuv/vmMH/pr9+/6DCeP/Hv6T/5qy//+Cotv/jl7L/+53R//yV0P/us+X/5cjz/+2v5P/8otb//ZTP//2i1f/0q9j/96LX//uSzv/5s9j/+6DS//mi0P+8uJP/17q0//6Ozv/9kc3+/KjXdQAAAAAAAAAAAAAAAAAAAAAAAAAAba/8DJvJ/aep0f//qdH//4nA//+Iv///g7z//366//+Hv///h7///325//+Iv///p9D//366//98t9//pLiD/9yttP/SuK3/w8ei/9u5uP/zl8X//oXK//6V0f/+lND//ZnR//2Uz//6ldH/4qbm/+Or5//4rt3//rLc//2p2P/xt97/8szs//yRzv/7sN3/9J3Z/+ms1P/OtK3/85fF//2Z0v/9odXE+5fPGQAAAAAAAAAAAAAAAAAA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//+Au///j8P//5HC4f+oyqv/oMav/7DEj/+6t5H/q8eG/67LjP+qyIb/z7ao/9q3tv++w5v/4qjA/9+q6f/grer/3qrp/+ae4P/0oNiX/4/BAwAAAAD///kDuNn/npLF//+Iv///f7r//4a+//+Dvf//jcL//3i3//92tv//fbr//4e///+Wx///icD9/4G7+P+cxsb/qsiE/7nNmv+mw3//osN6/6nGhP/Ds5r/yLmi/96ks//8gMX/6JfY/92r6//snNv/9ZXT//qX0ur9gsU2AAAAAOHv/iKZyP/ZkMP//5fH//+Kwf//g73//4vB//+Kwf//iL///3+7//+KwP//hL3//47D//+ax+f/ib7l/6XKx/+vy47/wMOd/6q/gv+iw3v/oMN4/6LFe/+kxX3/tbuN/9Kvqf/Ussb/56Pi//yY0f/treP/6Kjj//Gf2oQAAAAA0eb+VYvB//aJwP//j8P//3+7//+Wx///j8P//4W+//+JwP//fLn//3m3//97uP//hr7//5nF1f+gyMn/jL7V/7LOl/+ow4H/0q6p/+Cwuv+7vJT/s8SP/6/JjP+1wo//vMuc/7K6kv/xjsf//YPI//qT0P/rod//4K7qwei06A+ey/+Kj8P//4W9//+Au///gbz//4nA//+Fvv//hb7//366//9+uv//dbX//4W+//+Hv///ir/v/6PIs/+ZwJH/pcV9/7HIjf/asLP/+5XO/9Ouqv/VrKv/3re6/+imv//Yqa3/t72P/96otP/9k9D/8KHc//Si2v/4pNnm+JfSMZnI/rSMwf//hr7//3m3//9+uv//fLn//4a///+Buvv/obvW/4q78v9tsf//gLv//4S9//+Fvff/q8y0/6LEeP+kxH3/t9CY/6vIh//bsLT/8KvL//Odx//+odf/+arc//yY0P/zm8b/85PD//2c0//truT/4q3p/+uv5fnyqd5rp8/+2JTG//+YyP//kMP//4G8//9ztP//fLn//4/D/v+gwuf/grbx/3Cy/v90tf7/fbn//5LE8v+Wwb7/tbyQ/7HJjv+xzZD/wLqZ/+WZtv/+k9D//pLP//yx3f/3x+n//afY//2Pzv/9i8z/+o/J//aWzf/luOn/6bHm//qX0bGRxP/jm8b4/3q1+P+Evf//nMr//3a2//98uf//jMH//3+7//9usv//ib34/43A+v+Nwv3/rs67/7XOmP+/uJb/tMuS/7DPj//br7T/+IzG/+yfwf/zmMX//aPW//694//9otb//YjK//2NzP/psMT/7qXG//Wy2P/5o9X//aPVmoS+/+emyvL/irz0/4K8//+Uxv//frr//4W+//+Iv///dLT//3W2//+Mwv//k8X//5bG9/+407b/wced/9LCsf/avbj/0Lqs/8LHof+/wpv/tsKQ/82xpP/1msj/6rzK//Omy//9ldD//ZvT//mUy//8odT//prT//2Y0f/8qtmGjML/2pjI//+RxP//i8H//4S9//97uP//icD//4vB//+KwP//ksT//4/D//+dy///l8Xa/7TNl//TsKv/4bG8/9e1sv/vnML/usqZ/7DPj/+5vZL/5Zy3/92rtP/A0KL/ycur/+ymxP/2pM7/9p/M//qt1//socP/7KnV//ii2JiIv//FgLv//4e///+Nwv//fLn//5HE//+JwP//fbr//325//9/uv//frr//47D//+XxND/sMyM/6nIhf+vzI3/xr+j/+qrw//Quqz/xb+h/7/Inv/Dx6L/yMKx/9XEzP/UsK3/4KW1/8iynv+8wJf/wrOY/7q+lP/aucX896nfbm6x/3iFvv/9erj//3W1//+Dvf//ps///4vB//+Hv///ksX//326//+JwP//ib7k/6nLrf+/xZv/wrea/9ezsf/YubX/xsCi//OZxf/1mcj/8p3H/+qjv//wn83/86jg//yUz//9j83/8J3E/96stv/VtL7/3bvS//ac0Nf+mdEugbv+R43C//KEvf//dLX//3q4//+Kwf//g73//4vB//+JwP//mMj//5DD8/+jx6b/rMmG/9Kmpf/1m8n/86XL/++kxv/Nxq7/9JvI//SZx//9kc7//prT//WTzf/ypNP/96HY//Gh3P/9ltH//JTS/+aw6v/hu/L/6bTpfQAAAAB2tf0dnMr/1JXG//+SxP//hL3//4e///+Hv///fLn//3y5//+Mwf//j8Da/7POlP+vy43/u8eZ/7u+lf+9wJn/4LK7/8HCnf/vqMj/4qu6//ya0f/9j83/8ZrP/+uj0f/5ktL/7bTd//mh0f/1odr/4rju/+bB8P3jxPJbAAAAAG6u9gSRxP6jocz//4a+//98uf//gbv//3m4//92tv//hr7//4S9//+kzOj/tdCY/7e4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" />
    </head>

    <body>
    <h2 style="text-align: center;"> Cell communities clustering report</h2>
    <br>
    <table class="center">
        <thead>
            <tr>
                <th>Parameters used</th>
            </tr>
            <tr>
                <th><hr></th>
            </tr>
        </thead>
        <tbody>
        <tr>
            <td class="rightpad">Annotation</td>
            <td>{args.annotation}</td>
        </tr>
        <tr>
            <td class="rightpad">Resolution</td>
            <td>{args.resolution}</td>
        </tr>
        <tr>
            <td class="rightpad">Spot size</td>
            <td>{args.spot_size}</td>
        </tr>
        <tr>
            <td class="rightpad">Total cells normalization</td>
            <td>{args.total_cell_norm}</td>
        </tr>
        <tr>
            <td class="rightpad">Downsample rate</td>
            <td>{args.downsample_rate}</td>
        </tr>
        <tr>
            <td class="rightpad">Entropy threshold</td>
            <td>{args.entropy_thres}</td>
        </tr>
        <tr>
            <td class="rightpad">Scatteredness threshold</td>
            <td>{args.scatter_thres}</td>
        </tr>
        <tr>
            <td class="rightpad">Window sizes</td>
            <td>{args.win_sizes}</td>
        </tr>
        <tr>
            <td class="rightpad">Sliding steps</td>
            <td>{args.sliding_steps}</td>
        </tr>
        </tbody>
    </table>

    <h1>Slices: </h1>
    <hr>
    <div>
        {get_slices_content(args.out_path, args.plotting)}
    </div>
    <h1>Overall: </h1>
    <div>
        {get_overall_content(args.out_path, args.plotting)}
    </div>
    <p>Thank you for your attention. Place sth else here</p>
    </body>

    </html>
    '''

    with open(f"{args.out_path}/Report.html", "w") as f:
        f.write(htmlstr)
