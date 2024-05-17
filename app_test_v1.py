import io
import dash
from dash.dependencies import Input, Output, State
from dash import dcc
from dash import html
import pandas as pd
import base64
import plotly.express as px
import json


#username_password_pair=[['siddharth','evl123'],['username','password']]

app = dash.Dash()
#auth=app.dash_auth.BasicAuth(app, username_password_pair)
server=app.server

app.layout = html.Div(children=[
    html.Div(className='div-level-1', children=[
        html.H1(
            children='VENOME- A proteome analysis tool for venom characterization',
            className='header'
        ),
        html.Div(className='div-level-2', children=[
            html.Div(children=[
                html.Div(children=[
                    'Software used for proteome analysis',
                    dcc.Dropdown(options=['PEAKS version X', 'Proteome Discoverer', 'MaxQuant'],
                                 value='PEAKS version X',
                                 id='software-used'
                                 )]
                ),
                html.Div([
                    'Upload the csv or excel result file',
                    dcc.Upload(
                        id='upload-data',
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select Files', className='hover-text')
                        ]),
                        style={
                            'width': '90%',
                            'height': '60px',
                            'lineHeight': '60px',
                            'borderWidth': '1px',
                            'borderStyle': 'dashed',
                            'borderRadius': '5px',
                            'textAlign': 'center',
                            'margin': '5px 0px'
                        },
                        multiple=False
                    ),
                    html.Div(id='uploaded-file-confirmation')
                ], style={'padding': '10px 0px'}),
                html.Div([
                    dcc.Textarea(
                        id='input-peakarea',
                        placeholder='Enter the relative fraction data',
                        style={'height': '200px', 'padding': '1%', 'width': '90%'}
                    )
                ], style={'padding': '10px 0px'}
                ),
                html.Div(
                    children=[html.Button("Compute Relative Abundance", id="compute-button", className='button')],
                    className='form-control'
                ),
            ], className='div-level-3'),
            html.Div(className='div-level-4', children=[
                dcc.Graph(id='pie-chart'),
                html.Div(children=[
                    html.Button("Download CSV", id="download-button", className='button', style={'display': 'none'}),
                    dcc.Download(id="download-tox-abundance")
                ], style={'display': 'flex', 'alignItems': 'center', 'justifyContent': 'center'}),
                dcc.Store(id='tox-abundance-store', data=None)
            ])
        ])
    ])
])


@app.callback(
    Output('uploaded-file-confirmation', 'children'),
    Input('upload-data', 'filename')
)
def update_output_data(filename):
    if filename is not None:
        return '{} uploaded successfully'.format(filename)


# Helper function to read uploaded CSV file
def parse_contents(contents):
    content_type, content_string = contents.split(',')
    decoded = io.StringIO(base64.b64decode(content_string).decode('utf-8', errors='replace'))
    df = pd.read_csv(decoded)
    return df


def input_data_conversion(input_data):
    data_tuples = list(input_data.split('\n'))
    area_file = pd.DataFrame(data_tuples, columns=['proportions']).dropna().astype(float)
    return area_file


def compute_toxin_abundance(peak_file, area_file):
    # Take unique toxins only
    peak_file_unique = peak_file.drop_duplicates(subset='-10lgP')
    toxin_file = pd.read_csv("assets/tox-compiled.csv")
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
        proportion_peak_file[new_col1] = normalized_peak_file[i] * area_file.iloc[j - 1, 0]

    final_analysis_file = pd.concat([normalized_peak_file, proportion_peak_file.fillna(0)], axis=1)
    column_sum = proportion_peak_file.fillna(0).sum(axis=1)
    final_analysis_file['total_sum'] = column_sum
    # Compute toxin abundance
    tox_abundance = final_analysis_file.groupby('toxin_family')['total_sum'].sum()
    return tox_abundance


@app.callback(
    [Output('pie-chart', 'figure'),
     Output('download-button', 'style')],
    [Input('compute-button', 'n_clicks')],
    [State('upload-data', 'contents'),
     State('input-peakarea', 'value')]
)
def display_tox_abundance(n_clicks, peak_contents, area_contents):
    if n_clicks is None:
        return {}, {'display': 'none'}
    try:
        # Read uploaded CSV files
        peak_df = parse_contents(peak_contents)
        area_df = input_data_conversion(area_contents)
        # Call the compute_toxin_abundance function on button click
        tox_abundance = compute_toxin_abundance(peak_df, area_df)

        # Create pie chart
        fig = px.pie(tox_abundance, names=tox_abundance.index, values='total_sum', hole=.3,
                     title='Relative Abundance of toxin families')

        # Set labels only for slices greater than 1%
        fig.update_traces(textinfo='label+percent',
                          textposition='inside',
                          hoverinfo='label+percent+value')

        # Filter for values greater than 1% for labeling
        values = tox_abundance.values
        labels = [f'{label}: {value * 100:.2f}%' if value > 0.01 else '' for label, value in
                  zip(tox_abundance.index, values)]
        fig.update_traces(text=labels)

        return fig, {'display': 'block'}

    except Exception as e:
        return f"Error : {str(e)}", {}


@app.callback(
    Output('tox-abundance-store', 'data'),  # Store the tox_abundance data in the Store component
    [Input('compute-button', 'n_clicks')],
    [State('upload-data', 'contents'),
     State('input-peakarea', 'value')]
)
def store_tox_abundance(n_clicks, peak_contents, area_contents):
    if n_clicks is None:
        raise dash.exceptions.PreventUpdate
    try:
        peak_df = parse_contents(peak_contents)
        area_df = input_data_conversion(area_contents)
        tox_abundance = compute_toxin_abundance(peak_df, area_df)
        print(tox_abundance)
        return tox_abundance.reset_index().to_json(orient='split')
    except Exception as e:
        return f"Error: {str(e)}"


@app.callback(
    Output('download-tox-abundance', 'data'),
    Input('download-button', 'n_clicks'),
    State('tox-abundance-store', 'data'),
    prevent_initial_call=True
)
def download_csv(n_clicks, tox_abundance_json):
    print(tox_abundance_json)
    if n_clicks is None:
        return None
    try:
        # Ensure the JSON data is correctly parsed and convert
        tox_abundance = pd.read_json(io.StringIO(tox_abundance_json),orient='split')
        print(tox_abundance)
        print(type(tox_abundance))
        # Create a CSV string from the DataFrame
        csv_string = tox_abundance.to_csv(index=False)
        print(csv_string)
        # Trigger the download
        return dict(content=csv_string, filename="tox_abundance.csv")
    except Exception as e:
        return f"Error : {str(e)}", {}


if __name__ == '__main__':
    app.run_server(debug=True)
