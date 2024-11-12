import os
from bs4 import BeautifulSoup
import argparse
import base64


opt_parser = argparse.ArgumentParser(
    description="Generate an HTML report including all HTML, SVG, and PNG files in a directory."
)
opt_parser.add_argument(
    "-d",
    "--directory",
    dest="directory",
    type=str,
    help="Directory containing the files to include in the report",
)
opt_parser.add_argument(
    "-o", "--output_path", dest="output_path", type=str, help="Output HTML file"
)

options = opt_parser.parse_args()

directory = options.directory
output_path = options.output_path


def convert_image_to_base64(image_path):
    with open(image_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode("utf-8")
    return encoded_string


def read_svg_file(svg_path):
    with open(svg_path, "r") as svg_file:
        svg_content = svg_file.read()
    return svg_content


# Function to generate the HTML report
def generate_html_report(directory, output_file):
    # Create a BeautifulSoup object to construct the HTML document
    html_content = """
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Report ribosomal RNA analysis</title>
        </head>
        <body>
        </body>
        </html>
        """
    soup = BeautifulSoup(html_content,'html.parser')
    
    # Add some CSS styling
    style_tag = soup.new_tag('style')
    style_tag.string = """
    body {
        font-family: Arial, sans-serif;
        margin: 0;
        padding: 0;
        justify-content: center;
        align-items: center;
    }
    .centered-image {
                    display: block;
                    margin: 0 auto;
                    text-align: center;
                    vertical-align: middle;
                }
    .header {
        background-color: #3b5861;
        padding: 10px;
        text-align: center;
        color: white;
        position: relative;
    }
    .logo {
        display: inline-block;
        vertical-align: middle;
    }
    .menu {
        display: inline-block;
        vertical-align: middle;
    }
    .menu a {
        color: white;
        padding: 14px 20px;
        text-decoration: none;
        display: inline-block;
    }
    .menu a:hover {
        background-color: #1c1c24;
    }
    .section {
    border: 1px solid #ddd;
    padding: 20px;
    margin: 20px 0;
    border-radius: 4px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    """
    #margin-left: 50px;
    soup.head.append(style_tag)
    body = soup.body

    # Add a header with a banner, a logo, and a menu bar
    header_tag = soup.new_tag('div', **{'class': 'header'})
    

    # Add a logo
    #logo_tag = soup.new_tag('img', src='logo.png', **{'class': 'logo', 'alt': 'Logo', 'width': '50', 'height': '50'})
    #header_tag.append(logo_tag)

    #Add a menu bar
    menu_tag = soup.new_tag('div', **{'class': 'menu'})
    header_tag.append(menu_tag)

    # Add menu items
    coverage_tag = soup.new_tag('a', href='#coverage')
    coverage_tag.string = 'Coverage reference'
    menu_tag.append(coverage_tag)

    template_association_tag = soup.new_tag('a', href='#template_association')
    template_association_tag.string = 'Coverage ribosomal fragments'
    menu_tag.append(template_association_tag)
    
    cut_sites_tag = soup.new_tag('a', href='#cut_sites')
    cut_sites_tag.string = 'Cut sites'
    menu_tag.append(template_association_tag)

    intensity_matrix_tag = soup.new_tag('a', href='#intensity_matrix')
    intensity_matrix_tag.string = 'Intensity Matrix'
    menu_tag.append(intensity_matrix_tag)
    
    intensity_clustering_tag = soup.new_tag('a', href='#intensity_clustering')
    intensity_clustering_tag.string = 'Intensity Clustering'
    menu_tag.append(intensity_clustering_tag)
    
    hdbscan_clustering_tag = soup.new_tag('a', href='#hdbscan_clustering')
    hdbscan_clustering_tag.string = 'HDBSCAN Clustering'
    menu_tag.append(hdbscan_clustering_tag)
    
    modification_detection_tag = soup.new_tag('a', href='#modification_detection')
    modification_detection_tag.string = 'Relative modification abundance'
    menu_tag.append(modification_detection_tag)
    
    body.append(header_tag)
    
    


########################
#                      #
#  Reference coverage  #
#                      #
########################   
        
    div_tag = soup.new_tag('div', style="text-align: center;")
    file_path = os.path.join(
        directory, "coverage_plots/coverage_total_sample_absolute.png"
    )
    if os.path.exists(file_path):
        coverage_title_div = soup.new_tag('div', **{'class': 'section'})
        coverage_title_tag = soup.new_tag('h1', id= 'coverage')
        coverage_title_tag.string = 'Coverage ribosomal RNA'
        coverage_title_div.append(coverage_title_tag)   
        coverage_plots_text_tag = soup.new_tag('p')
        coverage_plots_text_tag.string = 'Absolute (left) and relative read coverage (right) of the ribosomal 45S template.\
            The plot show the read abundance in comparison to the alignment position on the reference.'
        coverage_title_div.append(coverage_plots_text_tag)
        body.append(coverage_title_div)
        
        image_base64 = convert_image_to_base64(file_path)
        with open(file_path, "r") as f:
            # Create an <img> tag for PNG files
            img_tag = soup.new_tag(
                "img",
                src=f"data:image/png;base64,{image_base64}",
                alt="Taillength per intensity based cluster",
                attrs={"class": "centered-img"},
                width=800,
                height=800
            )
            body.append(img_tag)
            img_tag.wrap(div_tag)

            
    file_path = os.path.join(
        directory, "coverage_plots/coverage_total_sample_relative.png"
    )
    if os.path.exists(file_path):
        image_base64 = convert_image_to_base64(file_path)
        with open(file_path, "r") as f:
            # Create an <img> tag for PNG files
            img_tag = soup.new_tag(
                "img",
                src=f"data:image/png;base64,{image_base64}",
                alt="Taillength per intensity based cluster",
                attrs={"class": "centered-img"},
                width=800,
                height=800
            )
            body.append(img_tag)
            img_tag.wrap(div_tag)
            body.append(soup.new_tag("br"))
    
########################
#                      #
#  Fragment coverage   #
#                      #
########################           
    
    div_tag2 = soup.new_tag('div', style="text-align: center;")
    file_path = os.path.join(
        directory, "coverage_plots/coverage_fragments_absolute.png"
    )
    if os.path.exists(file_path):
        fragment_coverage_title_div = soup.new_tag('div', **{'class': 'section'})
        fragment_coverage_title_tag = soup.new_tag('h1', id='template_association')
        fragment_coverage_title_tag.string = 'Coverage of ribosomal templates'
        fragment_coverage_title_div.append(fragment_coverage_title_tag)
        fragment_coverage_plots_text_tag = soup.new_tag('p')
        fragment_coverage_plots_text_tag.string = 'Relative and absolute abundance of reads associated to known ribosomal templates.\
            The used templates have been determined by third party research projects. The association of reads to literature based\
            templates is performed by determining the minimal overlap of each template-read pair.\
            The template showing maximal overlap between all possible template-read pairs is determined as associated reference.'
        fragment_coverage_title_div.append(fragment_coverage_plots_text_tag)
        body.append(fragment_coverage_title_div)
        
        image_base64 = convert_image_to_base64(file_path)
        with open(file_path, "r") as f:
            # Create an <img> tag for PNG files
            img_tag = soup.new_tag(
                "img",
                src=f"data:image/png;base64,{image_base64}",
                alt="Taillength per intensity based cluster",
                attrs={"class": "centered-img"},
                width=800,
                height=800
            )
            body.append(img_tag)
            img_tag.wrap(div_tag2)
            
            
    file_path = os.path.join(
        directory, "coverage_plots/coverage_fragments_relative.png"
    )
    if os.path.exists(file_path):
        image_base64 = convert_image_to_base64(file_path)
        with open(file_path, "r") as f:
            # Create an <img> tag for PNG files
            img_tag = soup.new_tag(
                "img",
                src=f"data:image/png;base64,{image_base64}",
                alt="Taillength per intensity based cluster",
                attrs={"class": "centered-img"},
                width=800,
                height=800
            )
            body.append(img_tag)
            img_tag.wrap(div_tag2)
            body.append(soup.new_tag("br"))
            
    div_tag3 = soup.new_tag('div', style="text-align: center;")
    file_path = os.path.join(
        directory, "coverage_plots/coverage_fragments_absolute_all.png"
    )
    if os.path.exists(file_path):
        image_base64 = convert_image_to_base64(file_path)
        with open(file_path, "r") as f:
            # Create an <img> tag for PNG files
            img_tag = soup.new_tag(
                "img",
                src=f"data:image/png;base64,{image_base64}",
                alt="Taillength per intensity based cluster",
                attrs={"class": "centered-img"},
                width=800,
                height=800
            )
            body.append(img_tag)
            img_tag.wrap(div_tag3)
            body.append(soup.new_tag("br"))
            
    
            
    file_path = os.path.join(
        directory, "polyA_template_based/polyA_tails_intermediates_template.html"
    )
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            fragment_taillength_title_div = soup.new_tag('div')
            fragment_taillength_title_tag = soup.new_tag('h1')
            fragment_taillength_title_tag.string = 'Taillenght of ribosomal templates'
            fragment_taillength_title_div.append(fragment_taillength_title_tag)
            template_based_analysis_tag = soup.new_tag('p')
            template_based_analysis_tag.string = 'Interactive relative abundance plots of template associated reads.\
                The three plots show the reference based cut sites on the 45S reference (top), \
                the minimal and maximal of cut sites (center) and \
                the mean start and end sites (bottom) of the associated reads.\
                The mean taillength and standard deviation is shown for each template.\
                The violinplot shows the taillength distributions for each template.\
                The association of reads to literature based\
                templates is performed by determining the minimal overlap of each template-read pair.\
                The template showing a maximal overlap between all possible template-read pairs is determined as reference for the given read.\
                PolyA taillength is determined with the polya tail prediction function integrated in dorado.'
            fragment_taillength_title_div.append(template_based_analysis_tag)
            body.append(fragment_taillength_title_div)
            
            # Read the content of the HTML file and append it to the body
            file_content = BeautifulSoup(f.read(), "html.parser")
            body.append(file_content)

    file_path = os.path.join(
        directory, "polyA_template_based/polyA_tails_intermediates_min_max.html"
    )
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            # Read the content of the HTML file and append it to the body
            file_content = BeautifulSoup(f.read(), "html.parser")
            body.append(file_content)
    
    
    

            
    file_path = os.path.join(
        directory, "polyA_template_based/polyA_tails_intermediates_mean.html"
    )
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            # Read the content of the HTML file and append it to the body
            file_content = BeautifulSoup(f.read(), "html.parser")
            body.append(file_content)

    
    
    
    file_path = os.path.join(
        directory, "polyA_template_based/violinplot_taillength_per_intermediate.png"
    )
    if os.path.exists(file_path):
        image_base64 = convert_image_to_base64(file_path)
        with open(file_path, "r") as f:
            # Create an <img> tag for PNG files
            img_tag = soup.new_tag(
                "img",
                src=f"data:image/png;base64,{image_base64}",
                alt="Taillength per intensity based cluster",
                attrs={"class": "centered-img"},
                width=1000,
                height=1000
            )
            div_tag = soup.new_tag('div', style="text-align: center;")
            body.append(img_tag)
            img_tag.wrap(div_tag)
            body.append(soup.new_tag("br"))

########################
#                      #
#       Cut sites      #
#                      #
########################

    
    
    file_path = os.path.join(
        directory, "cut_site_plots/cut_sites.html"
    )
        
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            cut_sites_title_div = soup.new_tag('div', **{'class': 'section'})
            cut_sites_title_tag = soup.new_tag('h1',id='cut_sites')
            cut_sites_title_tag.string = 'Cut sites'
            cut_sites_title_div.append(cut_sites_title_tag)
            cut_sites_text_tag = soup.new_tag('p')
            cut_sites_text_tag.string = 'Cut sites occurring in significant abundance in relation\
                to the gaussian distribution of cut sites of literature based templates.\
                Relative cut site abundance has been normalized by the absolute number of reads aligning to the 45S reference.\
                For intersecting templates the significance threshold has been determined by the combined gaussian of both fragments.'
            cut_sites_title_div.append(cut_sites_text_tag)
            body.append(cut_sites_title_div)
            body.append(soup.new_tag("br"))
        
            # Read the content of the HTML file and append it to the body
            file_content = BeautifulSoup(f.read(), "html.parser")
            body.append(file_content)
            body.append(soup.new_tag("br"))
            
            
########################
#                      #
#   Intensity Matrix   #
#                      #
########################      

    file_path = os.path.join(directory, "intensity_matrix/intensity_matrix.png")
    if os.path.exists(file_path):
        intensity_matrix_title_div = soup.new_tag('div', **{'class': 'section'})
        intensity_matrix_title_tag = soup.new_tag('h1',id='intensity_matrix')
        intensity_matrix_title_tag.string = 'Intensity matrix'
        intensity_matrix_title_div.append(intensity_matrix_title_tag)
        intensity_matrix_tag = soup.new_tag('p')
        intensity_matrix_tag.string = 'The intensity matrix plot shows the min-max\
            normalized amount of reads sharing the same start and end points on the 45S reference after alignment (top).\
            The interactive plot (bottom) represents a sub sampled version\
            of the intensity matrix showing the 15000 most abundant clusters.'
        intensity_matrix_title_div.append(intensity_matrix_tag)
        body.append(intensity_matrix_title_div)   
        
        image_base64 = convert_image_to_base64(file_path)
        with open(file_path, "r") as f:

            # Create an <img> tag for PNG files
            img_tag = soup.new_tag(
                "img",
                src=f"data:image/png;base64,{image_base64}",
                alt="Intensity Matrix of clusters sorted by start and endpoints of alignment",
                attrs={"class": "centered-img"},
                width=800,
                height=800
            )
            div_tag = soup.new_tag('div', style="text-align: center;")
            body.append(img_tag)
            img_tag.wrap(div_tag)
            body.append(soup.new_tag("br"))
    
    file_path = os.path.join(
        directory, "intensity_matrix/intensity_matrix.html"
    )
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            file_content = BeautifulSoup(f.read(), "html.parser")
            body.append(file_content)



########################
#                      #
# Intensity Clustering #
#                      #
########################   
    
    
    file_path = os.path.join(
        directory, "polyA_intensity_based_clusters/polyA_tails_clustering.html"
    )
    if os.path.exists(file_path):
        intensity_clustering_title_div = soup.new_tag('div', **{'class': 'section'})
        intensity_clustering_title_tag = soup.new_tag('h1',id='intensity_clustering')
        intensity_clustering_title_tag.string = 'Intensity Clustering'
        intensity_clustering_title_div.append(intensity_clustering_title_tag)
        intensity_clustering_tag = soup.new_tag('p')
        intensity_clustering_tag.string = 'Interactive plot showing the min max normalized abundance and position on 45S reference of read clusters\
            determined by reads sharing the same start and end sites on the 45S matrix.\
            The 300 biggest clusters are visualized.\
            Additionally, the mean and stdd of polyA-tails of the clusters are shown.\
            The taillength is determined by the polya tail prediction function integrated in dorado.\
            For validation of the clustering approach the 300 biggest clusters are shown as squares above the instensity matrix\
            ranging from start to end of reference on the coordinate system. The Upper left corner should allocate with a cluster showing a high intensity or read density.'
        intensity_clustering_title_div.append(intensity_clustering_tag)
        body.append(intensity_clustering_title_div)    
        
        with open(file_path, "r") as f:
            # Read the content of the HTML file and append it to the body
            file_content = BeautifulSoup(f.read(), "html.parser")
            body.append(file_content)
    
    file_path = os.path.join(
        directory, "fragment_analysis_intensity/intensity_matrix.png"
    )
    if not os.path.exists(file_path):
        file_path = os.path.join(
        directory, "fragment_intensity_analysis/intensity_matrix.png"
    )
    if os.path.exists(file_path):
        image_base64 = convert_image_to_base64(file_path)
        with open(file_path, "r") as f:

            # Create an <img> tag for PNG files
            img_tag = soup.new_tag(
                "img",
                src=f"data:image/png;base64,{image_base64}",
                alt="Clusters detected on intensity matrix",
                attrs={"class": "centered-img"},
                width=800,
                height=800
            )
            div_tag = soup.new_tag('div', style="text-align: center;")
            body.append(img_tag)
            img_tag.wrap(div_tag)

    file_path = os.path.join(
        directory,
        "polyA_intensity_based_clusters/violinplot_taillength_per_cluster.png",
    )
    if os.path.exists(file_path):
        image_base64 = convert_image_to_base64(file_path)
        with open(file_path, "r") as f:
            # Create an <img> tag for PNG files
            img_tag = soup.new_tag(
                "img",
                src=f"data:image/png;base64,{image_base64}",
                alt="Taillength per intensity based cluster",
                attrs={"class": "centered-img"},
                width=1000,
                height=1000
            )
            div_tag = soup.new_tag('div', style="text-align: center;")
            body.append(img_tag)
            img_tag.wrap(div_tag)
            body.append(soup.new_tag("br"))
            
    
########################
#                      #
#  HDBSCAN Clustering  #
#                      #
########################

    
    
    file_path = os.path.join(
        directory, "polyA_hdbscan_based_clusters/polyA_tails_clustering.html"
    )
    
    
    if os.path.exists(file_path):
        hdbscan_clustering_title_div = soup.new_tag('div', **{'class': 'section'})
        hdbscan_clustering_title_tag = soup.new_tag('h1', id='hdbscan_clustering')
        hdbscan_clustering_title_tag.string = 'HDBSCAN Clustering'
        hdbscan_clustering_title_div.append(hdbscan_clustering_title_tag)
        hdbscan_clustering_tag = soup.new_tag('p')
        hdbscan_clustering_tag.string = 'Interactive plot showing the min max normalized abundance\
            and position on 45S reference of read clusters\
            determined by density based clustering approach using the intensity matrix as template.\
            The 300 biggest clusters are visualized.\
            Additionally, the mean and stdd of polyA-tails of the clusters are shown.\
            The taillength is determined by the polya tail prediction function integrated in dorado.\
            For validation of the clustering approach the 300 biggest clusters are shown as squares above the instensity matrix\
            ranging from start to end of reference on the coordinate system. The Upper left corner should allocate with a cluster showing a high intensity or read density.'
        hdbscan_clustering_title_div.append(hdbscan_clustering_tag)
        body.append(hdbscan_clustering_title_div)

        with open(file_path, "r") as f:
            # Read the content of the HTML file and append it to the body
            file_content = BeautifulSoup(f.read(), "html.parser")
            body.append(file_content)

    file_path = os.path.join(
        directory, "fragment_analysis_hdbscan/intensity_matrix.png"
    )
    if not os.path.exists(file_path):
        file_path = os.path.join(
        directory, "fragment_hdbscan_analysis/intensity_matrix.png"
    )
    if os.path.exists(file_path):
        image_base64 = convert_image_to_base64(file_path)
        with open(file_path, "r") as f:

            # Create an <img> tag for PNG files
            img_tag = soup.new_tag(
                "img",
                src=f"data:image/png;base64,{image_base64}",
                alt="Clusters detected on intensity matrix",
                attrs={"class": "centered-img"},
                width=800,
                height=800
            )
            div_tag = soup.new_tag('div', style="text-align: center;")
            body.append(img_tag)
            img_tag.wrap(div_tag)


    file_path = os.path.join(
        directory, "polyA_hdbscan_based_clusters/violinplot_taillength_per_cluster.png"
    )
    if os.path.exists(file_path):
        image_base64 = convert_image_to_base64(file_path)
        with open(file_path, "r") as f:

            # Create an <img> tag for PNG files
            img_tag = soup.new_tag(
                "img",
                src=f"data:image/png;base64,{image_base64}",
                alt="Taillength per hdbscan based cluster",
                attrs={"class": "centered-img"},
                width=1000,
                height=1000
            )
            div_tag = soup.new_tag('div', style="text-align: center;")
            body.append(img_tag)
            img_tag.wrap(div_tag)
            body.append(soup.new_tag("br"))
        
    

########################
#                      #
#Modification detection#
#                      #
########################

    
    file_path = os.path.join(
        directory, "modification_plots/relative_pseU_modification_abundance.html"
    )
    if os.path.exists(file_path):
        modification_title_div = soup.new_tag('div', **{'class': 'section'})
        modification_title_tag = soup.new_tag('h1',id='modification_detection')
        modification_title_tag.string = 'Modifications on ribosomal template'
        modification_title_div.append(modification_title_tag)
        modification_text_tag = soup.new_tag('p')
        modification_text_tag.string = 'Relative modification abundance plots for pseU (top) and m6A (bottom) modifications over the 45S reference.\
        The modification detection is performed with the integrated dorado models for modification detection.\
        In the pseU plot the systematic C/U missmatch reported for pseU is added, since the dorado models are not taking missbasecalled bases in consideration.' 
        modification_title_div.append(modification_text_tag)
        body.append(modification_title_div)
        
        with open(file_path, "r") as f:
            # Read the content of the HTML file and append it to the body
            file_content = BeautifulSoup(f.read(), "html.parser")
            body.append(file_content)
    
    file_path = os.path.join(
        directory, "modification_plots/relative_m6A_modification_abundance.html"
    )
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            # Read the content of the HTML file and append it to the body
            file_content = BeautifulSoup(f.read(), "html.parser")
            body.append(file_content)
            
    
    with open(output_file, "w") as f:
        f.write(soup.prettify())
    
print(f"HTML report generated and saved to {output_path}rRNA_report.html")
generate_html_report(directory, f"{output_path}rRNA_report.html")
