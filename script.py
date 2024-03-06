import io
import dash
from dash.dependencies import Input, Output
from dash import dcc
from dash import html
import pandas as pd
import base64
import plotly.express as px

app = dash.Dash()
server = app.server

app.layout = html.Div([
    html.Div([
        'Upload the exported file',
        dcc.Upload(
            id='upload-data-peak',
            children=[
                html.Button('Upload File'),
                html.Div(id='uploaded-raw-file-label')
            ])
    ]),
    html.Div([
        'Upload the HPLC peak area file',
        dcc.Upload(
            id='upload-data-area',
            children=[
                html.Button('Upload HPLC File'),
                html.Div(id='uploaded-HPLC-file-label')
            ])
    ]),
    html.Div([
        'Upload the toxin annotation file',
        dcc.Upload(
            id='upload-data-toxin',
            children=[
                html.Button('Upload toxin annotation File'),
                html.Div(id='uploaded-toxin-file-label')
            ])
    ]),
    html.Br(),
    html.Button("Compute Relative Abundance", id="compute-button"),
    html.Br(),
    dcc.Graph(id='pie-chart')
])


@app.callback(
    Output('uploaded-raw-file-label', 'children'),
    Input('upload-data-peak', 'filename')
)
def update_output_data(filename):
    if filename is not None:
        return '{} uploaded successfully'.format(filename)


@app.callback(
    Output('uploaded-HPLC-file-label', 'children'),
    Input('upload-data-area', 'filename')
)
def update_output_data_hplc(filename):
    if filename is not None:
        return '{} uploaded successfully'.format(filename)


@app.callback(
    Output('uploaded-toxin-file-label', 'children'),
    Input('upload-data-toxin', 'filename')
)
def update_output_data_hplc(filename):
    if filename is not None:
        return '{} uploaded successfully'.format(filename)


# Helper function to read uploaded CSV file
def parse_contents(contents):
    content_type, content_string = contents.split(',')
    decoded = io.StringIO(base64.b64decode(content_string).decode('utf-8'))
    df = pd.read_csv(decoded)
    return df


def compute_toxin_abundance(peak_file, area_file, toxin_file):
    # Take unique toxins only
    peak_file_unique = peak_file.drop_duplicates(subset='-10lgP')

    # Merge dataframes
    merge_file = pd.merge(peak_file_unique, toxin_file, how='inner', on='Accession')

    # Filter columns
    area_columns = merge_file.filter(like='Area')
    selected_columns = ['Accession', 'Description', 'toxin_family'] + list(area_columns.columns)
    filtered_peak_file = merge_file[selected_columns].fillna(0)

    # Adding normalized data
    normalized_peak_file1 = pd.DataFrame()
    for i in area_columns.columns:
        normalized_area_columns = f"{i}_N"
        normalized_peak_file1[normalized_area_columns] = area_columns[i] / area_columns[i].sum()

    normalized_peak_file = pd.concat([filtered_peak_file, normalized_peak_file1], axis=1).fillna(0)

    # Adding proportion column
    proportion_peak_file = pd.DataFrame()
    j = 0
    for i in normalized_peak_file1.columns:
        j += 1
        new_col1 = f"{i}_P"
        proportion_peak_file[new_col1] = normalized_peak_file1[i] * area_file.iloc[j - 1, 0] * area_file.iloc[j - 1, 1]

    final_analysis_file = pd.concat([normalized_peak_file, proportion_peak_file.fillna(0)], axis=1)
    column_sum = proportion_peak_file.fillna(0).sum(axis=1)
    final_analysis_file['total_sum'] = column_sum

    # Compute toxin abundance
    tox_abundance = final_analysis_file.groupby('toxin_family')['total_sum'].sum()
    return tox_abundance


@app.callback(
    Output('pie-chart', 'figure'),
    [Input('compute-button', 'n_clicks')],
    [Input('upload-data-peak', 'contents'),
     Input('upload-data-area', 'contents'),
     Input('upload-data-toxin', 'contents')]
)
def display_tox_abundance(n_clicks, peak_contents, area_contents, toxin_contents):
    if n_clicks is None:
        raise dash.exceptions.PreventUpdate
    try:
        # Read uploaded CSV files
        peak_df = parse_contents(peak_contents)
        area_df = parse_contents(area_contents)
        toxin_df = parse_contents(toxin_contents)
        # Call the compute_toxin_abundance function on button click
        tox_abundance = compute_toxin_abundance(peak_df, area_df, toxin_df)
        # Create pie chart
        fig = px.pie(tox_abundance, names=tox_abundance.index, values='total_sum',
                     title='Relative Abundance')
        return fig
    except Exception as e:
        return f"Error: {str(e)}", {}


if __name__ == '__main__':
    app.run_server(debug=True)

# %%
