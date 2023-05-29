import base64
import os
import subprocess

from collections import defaultdict
from .utils import timeit


BOXPLT_C_INDEX = 9
COLORPLT_C_INDEX = 15
CMIXT_C_INDEX = -5

def get_base64_encoded_img(path):
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
        footer {
            position: absolute;
            bottom: 0;
            width: 100%;
            height: 2.5rem;
            background-color:#FAEFFF;
            text-align: center;
            padding-bottom:20px
        }
        .content-wrap {
            padding-bottom: 2.5rem;
        }
        .page-container {
        position: relative;
        min-height: 100vh;
        }
        .pad20{
            padding: 0px 20px 0px;
        }
        .centered{
            width: 60%;
            margin: auto;
            text-align: center;
        }
        .centeredSmall{
            width: 70%;
            margin: auto;
            text-align: center;
        }
        .myButton {
            margin-top: 10px;
            margin-bottom:30px;
            box-shadow:inset 0px 1px 0px 0px #efdcfb;
            background-color:#5F0085;
            border-radius:6px;
            border:1px solid #c584f3;
            display:inline-block;
            cursor:pointer;
            color:#ffffff;
            font-family:Arial;
            font-size:15px;
            font-weight:bold;
            padding:6px 24px;
            text-decoration:none;
            text-shadow:0px 1px 0px #9752cc;
        }
        .myButton:hover {
            background:linear-gradient(to bottom, #bc80ea 5%, #dfbdfa 100%);
            background-color:#bc80ea;
        }
        .myButton:active {
            position:relative;
            top:1px;
        }
    '''

def all_slices_get_figure(project_root, figure_name, plotting_level):
    if plotting_level > 0:
        return get_base64_encoded_img(f"{project_root}/{figure_name}") if os.path.isfile(f"{project_root}/{figure_name}") else ""

def per_slice_content(path, plotting_level):
    content = ""
    for root, dirs, files in os.walk(path): 
        for name in dirs:
            content += get_table_plots(os.path.join(root, name))
    return content

def get_table_plots(path):
    for root, dirs, files in os.walk(path):
        for name in files:
            if name.startswith("celltype_table"):
                celltype_table = os.path.join(root, name)
            if name.startswith("cell_mixture_table"):
                cell_mixtures = os.path.join(root, name)
            if name.startswith("window_cell_num_hist"):
                hist_cell_number = os.path.join(root, name)
    return f'''
        <table>
            <thead>
                <tr>
                    <th colspan="3"><h3>slice: {os.path.basename(path)}</h3></th>
                </tr>
            </thead>
            <tbody>
            <tr>
                <td class="pad20" style="width:40%"><img src={get_base64_encoded_img(celltype_table)}></td>
                <td class="pad20"><img src={get_base64_encoded_img(cell_mixtures)}></td>
                <td class="pad20" style="width:20%"><img src={get_base64_encoded_img(hist_cell_number)}></td>
            </tr>
            </tbody>
        </table>
        <hr>
    '''

def per_community_content(path, plotting_level):
    content = ""
    cmixtures_dict = defaultdict(list)
    boxplots_dict = defaultdict(list)
    colorplot_dict = defaultdict(list)

    for root, dirs, files in os.walk(path): 
        for name in dirs:
            slice_path = os.path.join(path, name)
            for file in os.listdir(slice_path):
                if file.startswith("boxplot"):
                    boxplots_dict[file[BOXPLT_C_INDEX]].append(os.path.join(slice_path, file))
                if file.startswith("cmixtures"):
                    cmixtures_dict[file[CMIXT_C_INDEX]].append(os.path.join(slice_path, file))
                if file.startswith("colorplot"):
                    colorplot_dict[file[COLORPLT_C_INDEX]].append(os.path.join(slice_path, file))
    
    content += make_table(cmixtures_dict, 2, "Cell types that are present in each community") #TODO if needed additional comments what these plots represent
    content += make_table(boxplots_dict, 2, "Boxplots of cell types that are present in each community")
    content += make_table(colorplot_dict, 2, "RGB Colorplots")
    return content

def make_table(plot_dict, columns, comment):
    content = ""
    for _, plots in plot_dict.items():
        rows = ""
        plots = [f'<td class="shrink"><img src={get_base64_encoded_img(plot)}></td>' for plot in plots]
        for i in range(0,len(plots),columns):
            row = "".join(plots[i:i+columns])
            rows += f'<tr>{row}</tr>'
        content += f'''
        <table style="border: 1px solid black;">
            <tbody class="twobytwo">
            {rows}
            </tbody>
        </table>
        <hr>
        '''

    return f'''
        <table>
            <thead>
                <tr>
                    <td colspan="{columns}"><h4>{comment}</h4></td>
                </tr>
            </thead>
        </table>
        {content}
        <hr>
    '''

@timeit
def generate_report(args):
    if args.plotting == 0:
        return

    command = "python main.py "
    for k,v in vars(args).items():
        if v == None or k == 'out_path' or k == 'project_name' or k == "skip_stats" or k == "save_adata":
            continue
        if k == 'out_path_orig' or k == 'project_name_orig':
            k = k[:-5]
        command += f"--{k} {v} "
    commit_date = subprocess.check_output(r'git log --pretty=format:"%h%x09%x09%ad%x09%s" -n 1', shell=True).expandtabs()

    htmlstr = f'''
    <!DOCTYPE html>

    <html lang=”en”>
    <head>
    <style>
    {get_css()}
    </style>

    <script>
    let text = "{command}"
    const copyCommand = async () => {{
        try {{
        await navigator.clipboard.writeText(text);
        console.log('Command copied to clipboard');
        }} catch (err) {{
        console.error('Failed to copy: ', err);
        }}
    }}
    </script>

    <title>Cell Communities Report</title>
    <link rel="icon" type="image/png" sizes="16x16" href={get_base64_encoded_img(f"{os.path.dirname(__file__)}/assets/favicon.ico")}/>
    </head>

    <body>
    <div class="page-container">
        <div class="content-wrap">
            <div style="background-color:#FAEFFF;">
                <h1 style="text-align: center; margin:0"> Cell communities clustering report</h1>
                <br>
                <table class="center">
                    <thead>
                        <tr>
                            <th colspan="2"><h3>Parameters used</h3></th>
                        </tr>
                        <tr>
                            <th colspan="2"><hr></th>
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
                    <tr>
                    <td colspan="2"><button title="Command that is used to generate results is copied to clipboard" class="rightpad myButton" type="button" onclick="copyCommand()">
                    Copy command</button></td>
                    </tr>
                    </tbody>
                </table>
            </div>
            <br><hr>
            <div>
                <h2 style="text-align: center;">Overall</h2>
            </div>
            <div class="centered">
                <h3>Cell type annotation:</h3>
            </div>
            <div class="centered">
                <img title="Annotation column in .obs of anndata object" src={all_slices_get_figure(args.out_path, "cell_type_per_slice.png", args.plotting)}>
            </div>
            <div class="centered">
                <h3>Communities obtained:</h3>
            </div>
            <div class="centered">
                <img title="Result of CCD algorithm contained in .obs of anndata object" src={all_slices_get_figure(args.out_path, "clustering_per_slice.png", args.plotting)}>
            </div>
            <div class="centered">
                <h3>Cell mixtures for all slices:</h3>
            </div>
            <div class="centered">
                <img title="Showing the percentage of each cell type within obtained community and corresponding sums" src={all_slices_get_figure(args.out_path, "total_cell_mixtures_table.png", args.plotting)}>
            </div>
            <div class="centered">
                <h3>Cell type abundance:</h3>
            </div>
            <div class="centered">
                <img title="Percentage of specific cell types per slice" src={all_slices_get_figure(args.out_path, "cell_abundance_all_slices.png", args.plotting)}>
            </div>
            <div class="centered">
                <h3>Communities abundance:</h3>
            </div>
            <div class="centered">
                <img title="Percentage of cells in each community per slice" src={all_slices_get_figure(args.out_path, "cluster_abundance_all_slices.png", args.plotting)}>
            </div>
            <div class="centered">
                <h3>Cell percentage in communities:</h3>
            </div class="centered">
            <div class="centered">
                <img title="Percentage of cells in each community per slice" src={all_slices_get_figure(args.out_path, "cell_perc_in_community_per_slice.png", args.plotting)}>
            </div>
            <br><hr>
            <div>
                <h2 style="text-align: center;">Per slice</h2>
            </div>
            <div>
                {per_slice_content(args.out_path, args.plotting)}
            </div>
            <br><hr>
            <div>
                <h2 style="text-align: center;">Per community</h2>
            </div>
            <div class="centeredSmall">
                {per_community_content(args.out_path, args.plotting)}
            </div>
        </div>
        <footer><h4>Report created from commit: {commit_date}</h4></footer>
    </div>
    </body>
    </html>
    '''

    with open(f"{args.out_path}/report.html", "w") as f:
        f.write(htmlstr)
    
if __name__ == '__main__':
    #change .utils to utils
    import pickle 
    print(os.getcwd())
    file = open("core/args", "rb")
    args = pickle.load(file)
    file.close()
    generate_report(args)