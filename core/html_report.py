import base64
import os

from .utils import timeit


def get_base64_encoded_img(path):
    data_uri = base64.b64encode(open(path, 'rb').read()).decode('utf-8')
    return '<img src="data:image/png;base64,{0}">'.format(data_uri)

def load_favicon():
    path = f"{os.path.dirname(__file__)}/assets/favicon.ico"
    data_uri = base64.b64encode(open(path, 'rb').read()).decode('utf-8')
    return '"data:image/png;base64,{0}"'.format(data_uri)

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
    
    commit_date = os.system(r'git log --pretty=format:"%h%x09%x09%ad%x09%s" -n 1')

    htmlstr = f'''
    <!DOCTYPE html>

    <html lang=”en”>

    <head>
    <style>
    {get_css()}
    </style>

    <title>Cell Communities Report</title>
    <link rel="icon" type="image/png" sizes="16x16" href={load_favicon()}/>
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
    <p>{commit_date}</p>
    </body>

    </html>
    '''

    with open(f"{args.out_path}/report.html", "w") as f:
        f.write(htmlstr)
