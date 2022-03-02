import dash
from dash import dcc
from dash import html
from dash import Output, Input
import plotly.express as px
import dash_bootstrap_components as dbc
from dash import dash_table
from .plot_funcs import plot_data_seqdiv, plot_data_tajimas, plot_data_fst, plot_data_watt_thet
import pandas as pd
import allel
import numpy as np
import h5py
import sys
from .stats_funcs import pop_dict


def get_data(uid, r=4):
    # Get main data
    import_fields = list(pd.read_csv(f'website/uploads/{uid}_import_fields.csv',
                                     header=None)[0])

    main_data = allel.vcf_to_dataframe(f'website/uploads/{uid}_variants.vcf',
                                       fields=import_fields,
                                       alt_number=1)

    # Drop any gene columns if all data is NaN (i.e. no overlapping gene present in region):
    gene_cols = ['GENE1', 'GENE2']
    for g in range(len(gene_cols)):
        if (main_data.loc[:, gene_cols[g]].astype(str) == 'nan').all():
            main_data = main_data.drop(columns=[gene_cols[g]])

    # Rename GENE1 or GENE2 to GENE if only one remains in the query range:
    if 'GENE1' not in main_data.columns:
        main_data = main_data.rename(columns={"GENE2": "GENE"})
    elif 'GENE2' not in main_data.columns:
        main_data = main_data.rename(columns={"GENE1": "GENE"})

    # Rounding Data - Selecting fields based on import fields as this is static:
    main_data[import_fields[8::]] = main_data[import_fields[8::]].astype(
        str).astype(float).round(r)

    # Get summary CSVs
    summary_stats = pd.read_csv(f'website/uploads/{uid}_stats_df.csv').round(r)
    try:
        summary_fst = pd.read_csv(f'website/uploads/{uid}_fst_df.csv').round(r)
    except:
        summary_fst = None

    # Get seg_pos, ac_seg for plotting
    seg_pos = np.load(f'website/uploads/{uid}_seg_pos.npy')
    ac_seg = allel.AlleleCountsChunkedTable(
                                (h5py.File(f'website/uploads/{uid}_ac_seg.h5', 'r')
                                 )['ac_seg'])

    return main_data, summary_fst, summary_stats, seg_pos, ac_seg


def init_dashboard(server):
    """Create a Plotly Dash dashboard."""
    dash_app = dash.Dash(
        server=server,
        routes_pathname_prefix="/results/",
        external_stylesheets=[
            "/static/style.css",  # TODO: Figure out how to make style like rest of App
            "https://fonts.googleapis.com/css?family=Open+Sans:300,300i,400,400i,600,600i,700,700i|Raleway:300,300i,400,400i,500,500i,600,600i,700,700i|Poppins:300,300i,400,400i,500,500i,600,600i,700,700i",
            dbc.themes.BOOTSTRAP
        ],
    )

    dash_app.index_string = '''
            <!DOCTYPE html>
            <html>
                <head>
                    {%metas%}
                    <title>Query Results</title>
                    <link href="/static/img/qmul_64.png" rel="shortcut icon">
                    <link href="/static/img/qmul_icon.png" rel="apple-touch-icon">
                    {%css%}
                </head>
                <body>
                    {%app_entry%}
                    <footer>
                        {%config%}
                        {%scripts%}
                        {%renderer%}
                    </footer>
                </body>
            </html>
            '''

    dash_app.layout = html.Div([
        dcc.Location(id='url', refresh=False),
        html.Div(id='page-content')
    ])

    init_callbacks(dash_app)

    return dash_app.server


def display_data(uid):
    # Read in the data
    _, summary_fst, summary_stats, seg_pos, _ = get_data(uid)

    # Get the names of the populations the user has specified:
    user_pops = summary_stats['Population'].to_list()
    list_of_lists = [[key for key, value in pop_dict.items() if value == x] for x in user_pops]
    user_cond_list = [item for sublist in list_of_lists for item in sublist]

    # Create the App layout
    main_table = dbc.Card(
        [
            dbc.CardBody(
                [
                    html.Div(
                        [
                            html.Br(),
                            html.H5("Query Results", className="card-title",
                                    style={'font-weight': 'bold', 'display': 'inline-block'}),
                            dcc.Link(html.Button("Submit Another Query", className="query",
                                                 style={'display': 'inline-block'}),
                                     href='/server', refresh=True,
                                     ),

                        ],
                    )
                ]
            ),

            dbc.CardHeader(
                dbc.Tabs(
                    [
                        dbc.Tab(label="Allele Freq.", tab_id="allele_df"),
                        dbc.Tab(label="Genotype Freq.", tab_id="geno_df"),
                        dbc.Tab(label="Derived Allele Freq.", tab_id="daf_df")
                    ],
                    id="tabs",
                    active_tab="allele_df"
                )
            ),
            dbc.CardBody(html.Div(id="tab-content", className="p-4")),
        ], style={'border-color':'white'}
    )

    # Don't create fst table if user hasn't requested fst values:
    if summary_fst is not None:
        summary = dbc.Card(
            dbc.CardBody(
                [
                    html.Div([
                         html.H5("Summary Statistics", className="card-title",
                                 style={'font-weight': 'bold'}),
                         html.P(
                            "Summary statistics have been calculated over the entire range of the query."),
                         dbc.Col(create_small_table(summary_stats, 'stats')),
                        ], style={'display': 'inline-block', 'margin-right': '8em', "width": "44%"}),
                    html.Div([
                         html.H5("Fst - Population Comparisons", className="card-title",
                                 style={'font-weight': 'bold'}),
                         html.P(
                            "Fst values are shown for every possible population combination in the query."),
                         dbc.Col(create_small_table(summary_fst, "fst")),
                        ], style={'display': 'inline-block', "width": "44%"}),
                 ]
            ),
            color="white", outline=True,
        )
    else:
        summary = dbc.Card(
            dbc.CardBody(
                [
                    html.H5("Summary Statistics", className="card-title",
                            style={'font-weight': 'bold'}),
                    dbc.Col(create_small_table(summary_stats, 'stats')),
                ]
            ),
            color="white", outline=True
        )

    # Don't display charts if there are less than 10 variants:
    if len(seg_pos) > 10:
        charts = dbc.Card(
            dbc.CardBody(
                [
                    html.Br(),
                    html.H5("Summary Statistics Charts", className="card-title", style={'font-weight': 'bold'}),
                    html.P(
                        "Please choose an window size for the statistics to be averaged across and a step function to \
                         determine how far each window shifts (defaults set to 1,000bp and 100bp respectively). Please \
                         ensure both are lower than the number of SNPs in the data."),
                    dbc.Col(
                        [html.Div([
                            html.Div(
                                [
                                    html.H6("""Window Size (bp)""",
                                            style={'margin-right': '2em'})
                                ],
                            ),
                            dcc.Input(
                                id='window',
                                type='number',
                                value=1000,
                                min=1,
                                className="form-control w-auto",
                                style={"width": "50%", 'display': 'inline-block',
                                       'verticalAlign': "left", 'margin-right': '12em'}),
                            html.Div(
                                [
                                    html.H6("""Step (bp)""",
                                            style={'margin-right': '2em'})
                                ],
                            ),
                            dcc.Input(
                                id='step',
                                type='number',
                                value=100,
                                min=1,
                                className="form-control w-auto",
                                style={"width": "50%", 'display': 'inline-block', 'verticalAlign': "right"})
                        ],
                            style={'display': 'flex'}, ),
                            html.Br(),
                            html.Div([
                                html.Div(
                                    [
                                        html.H6("""First Population""",
                                                style={'margin-right': '2em'})
                                    ],
                                ),
                                dcc.Dropdown(user_cond_list,
                                             placeholder='-',
                                             id='pop_1',
                                             style={"width": "50%", 'display': 'inline-block', 'verticalAlign': "left"}),
                                html.Div(
                                    [
                                        html.H6("""Second Population""",
                                                style={'margin-right': '2em'})
                                    ],
                                ),
                                dcc.Dropdown(user_cond_list,
                                             placeholder='-',
                                             id='pop_2',
                                             style={"width": "50%", 'display': 'inline-block', 'verticalAlign': "right"})
                            ],
                                style={'display': 'flex'}, ),
                        ]),
                    html.Br(),
                    html.Div(id='graph_container',
                         children=[
                            dbc.Row(
                                [
                                    dbc.Col(dcc.Graph(id='seq_d_graph'), style={"width": "30%", 'margin-right': '6em'}),
                                    dbc.Col(dcc.Graph(id='taj_d_graph'), style={"width": "30%"}),
                                ]
                            ),
                            dbc.Row(
                                [
                                    dbc.Col(dcc.Graph(id='watt_thet_graph'), style={"width": "30%", 'margin-right': '6em'}),
                                    dbc.Col(dcc.Graph(id='fst_graph'), style={"width": "30%"}),
                                ],
                            )]
                             )
                ]
            )
            , color="white", outline=True
        )

    else:
        charts = dbc.Card(
            dbc.CardBody(
                [
                    html.Br(),
                    html.H5("Summary Statistics Charts", className="card-title", style={'font-weight': 'bold'}),
                    html.P(
                        "Summary statistics charts are not displayed for queries containing less than \
                        ten SNPs in the selected region.")
                    ]
            )
            , color="white", outline=True
        )

    return dbc.Container(
        [
            dbc.Row(
                [main_table]
            ),
            dbc.Row(
                [summary]
            ),
            dbc.Row(
                [charts]
            )
        ]
    )


def uid_from_url(url):
    """ Returning the file locations from the web_server.py redirect """
    file_loc = url.split("/")[-2:]
    sep = '/'
    uid = sep.join(file_loc)
    return uid


def init_callbacks(dash_app):
    @dash_app.callback(
        Output("page-content", "children"),
        Input("url", "pathname"))
    def render(url):
        uid = uid_from_url(url)
        return display_data(uid)

    @dash_app.callback(
        Output("graph_container", "style"),
        [Input('pop_1', 'value')],)
    def hide_graphs(p1):
        if p1 is not None:
            return {'display': 'block'}
        else:
            return {'visibility': 'hidden'}

    @dash_app.callback(
        Output("tab-content", "children"),
        [Input("url", "pathname"),
         Input("tabs", "active_tab")], )
    def update_tab_content(url, active_tab):

        uid = uid_from_url(url)
        df, _, _, _, _ = get_data(uid)

        # Listing the possible column values for each tab:
        a_cols = ['CHROM', 'POS', 'REF', 'ALT', 'GENE', 'GENE1', 'GENE2', 'RS_VAL', 'AF_AFR',
                  'AF_AMR', 'AF_EAS', 'AF_EUR', 'AF_SAS']
        g_cols = ['CHROM', 'POS', 'GF_HET_AFR', 'GF_HOM_REF_AFR', 'GF_HOM_ALT_AFR', 'GF_HET_AMR',
                  'GF_HOM_REF_AMR', 'GF_HOM_ALT_AMR', 'GF_HET_EAS', 'GF_HOM_REF_EAS',
                  'GF_HOM_ALT_EAS', 'GF_HET_EUR', 'GF_HOM_REF_EUR', 'GF_HOM_ALT_EUR',
                  'GF_HET_SAS', 'GF_HOM_REF_SAS', 'GF_HOM_ALT_SAS']
        d_cols = ['CHROM', 'POS', 'REF', 'ALT', 'AA', 'RS_VAL', 'DAF_AFR', 'DAF_AMR',
                  'DAF_EAS', 'DAF_EUR', 'DAF_SAS']

        # Returning the actual column values given by the user selection:
        df_col_names = df.columns
        a_freq_cols = [x for x in a_cols if x in df_col_names]
        g_freq_cols = [x for x in g_cols if x in df_col_names]
        daf_cols = [x for x in d_cols if x in df_col_names]

        if active_tab == 'allele_df':
            return dbc.Col(create_data_table(df[a_freq_cols]))
        if active_tab == 'geno_df':
            return dbc.Col(create_data_table(df[g_freq_cols]))
        if active_tab == 'daf_df':
            return dbc.Col(create_data_table(df[daf_cols]))
        return

    @dash_app.callback(
        Output('seq_d_graph', 'figure'),
        [Input("url", "pathname"),
         Input('window', 'value'),
         Input('step', 'value'),
         Input('pop_1', 'value'),
         Input('pop_2', 'value')],
    )
    def update_seq_graph(url, w, s, p1, p2):
        uid = uid_from_url(url)
        _, _, _, sp_df, ac_df = get_data(uid)

        if w and p1 is not None:
            seq_plot_df = plot_data_seqdiv(sp_df, ac_df, p1, p2, w, s)

            try:
                fig = px.line(seq_plot_df, x="chrom_pos", y="seq_div",
                              color="population", hover_name="population",
                              title='Nucleotide Diversity',
                              labels=dict(chrom_pos="Chromosome Position (bp)",
                                          seq_div="Nucleotide Diversity",
                                          population="Population"
                                          )
                              )
            except ValueError:
                fig = px.line(seq_plot_df, x="chrom_pos", y="seq_div",
                              color="population", hover_name="population",
                              title='Nucleotide Diversity',
                              labels=dict(chrom_pos="Chromosome Position (bp)",
                                          seq_div="Nucleotide Diversity",
                                          population="Population"
                                          )
                              )

            plotly_fmt(fig)
            return fig
        else:
            return "Please select window size and populations to compare"

    @dash_app.callback(
        Output('watt_thet_graph', 'figure'),
        [Input("url", "pathname"),
         Input('window', 'value'),
         Input('step', 'value'),
         Input('pop_1', 'value'),
         Input('pop_2', 'value')],
    )
    def update_wt_graph(url, w, s, p1, p2):
        uid = uid_from_url(url)
        _, _, _, sp_df, ac_df = get_data(uid)

        if w and p1 is not None:
            seq_plot_df = plot_data_watt_thet(sp_df, ac_df, p1, p2, w, s)

            try:
                fig = px.line(seq_plot_df, x="chrom_pos", y="watt_thet",
                              color="population", hover_name="population",
                              title='Watterson\'s Theta',
                              labels=dict(chrom_pos="Chromosome Position (bp)",
                                          watt_thet="Watterson's Theta",
                                          population="Population"
                                          )
                              )
            except ValueError:
                fig = px.line(seq_plot_df, x="chrom_pos", y="watt_thet",
                              color="population", hover_name="population",
                              title='Watterson\'s Theta',
                              labels=dict(chrom_pos="Chromosome Position (bp)",
                                          watt_thet="Watterson's Theta",
                                          population="Population"
                                          )
                              )

            plotly_fmt(fig)
            return fig
        else:
            return "Please select window size and populations to compare"

    @dash_app.callback(
        Output('taj_d_graph', 'figure'),
        [Input("url", "pathname"),
         Input('window', 'value'),
         Input('step', 'value'),
         Input('pop_1', 'value'),
         Input('pop_2', 'value')]
    )
    def update_taj_graph(url, w, s, p1, p2):
        uid = uid_from_url(url)
        _, _, _, sp_df, ac_df = get_data(uid)

        if w and p1 is not None:
            taj_plot_df = plot_data_tajimas(sp_df, ac_df, p1, p2, w, s)

            try:
                fig = px.line(taj_plot_df, x="chrom_pos", y="tajima_d",
                              color="population", hover_name="population",
                              title='Tajima\'s D',
                              labels=dict(chrom_pos="Chromosome Position (bp)",
                                          tajima_d="Tajima's D",
                                          population="Population"
                                          )
                              )
            except ValueError:
                fig = px.line(taj_plot_df, x="chrom_pos", y="tajima_d",
                              color="population", hover_name="population",
                              title='Tajima\'s D',
                              labels=dict(chrom_pos="Chromosome Position (bp)",
                                          tajima_d="Tajima's D",
                                          population="Population"
                                          )
                              )

            plotly_fmt(fig)
            return fig
        else:
            return "Please select window size and populations to compare"

    @dash_app.callback(
        Output('fst_graph', 'figure'),
        [Input("url", "pathname"),
         Input('window', 'value'),
         Input('step', 'value'),
         Input('pop_1', 'value'),
         Input('pop_2', 'value')]
    )
    def update_fst_graph(url, w, s, p1, p2):
        uid = uid_from_url(url)
        _, fst_df, _, sp_df, ac_df = get_data(uid)

        # If user hasn't selected Fst then fst_df will be None:
        if fst_df is not None:
            if w and p1 and p2 is not None:
                fst_plot_df = plot_data_fst(sp_df, ac_df, p1, p2, w, s)

                try:
                    fig = px.line(fst_plot_df, x="chrom_pos", y="fst",
                                  color="population", hover_name="population",
                                  color_discrete_sequence=["darkcyan"],
                                  title='Fst (Hudson\'s)',
                                  labels=dict(chrom_pos="Chromosome Position (bp)",
                                              fst="Fst",
                                              population="Populations"
                                              )
                                  )
                except ValueError:
                    fig = px.line(fst_plot_df, x="chrom_pos", y="fst",
                                  color="population", hover_name="population",
                                  color_discrete_sequence=["darkcyan"],
                                  title='Fst (Hudson\'s)',
                                  labels=dict(chrom_pos="Chromosome Position (bp)",
                                  fst="Fst",
                                  population="Populations"
                                  )
                                  )

                plotly_fmt(fig)
                return fig
            else:
                return "Please select window size and populations to compare"
        else:
            return
    return


def plotly_fmt(fig):
    return fig.update_layout(
        xaxis=dict(
            showline=True,
            showgrid=False,
            showticklabels=True,
            linecolor='rgb(204, 204, 204)',
            linewidth=2,
            ticks='outside',
            tickfont=dict(
                family='sans-serif',
                size=12,
                color='rgb(82, 82, 82)',
            )
        ), yaxis=dict(
            showgrid=False,
            showline=True,
            showticklabels=True,
            linecolor='rgb(204, 204, 204)',
            linewidth=2,
            ticks='outside',
            tickfont=dict(
                family='sans-serif',
                size=12,
                color='rgb(82, 82, 82)',
            ),
        ),
        width=650,
        height=300,
        margin=dict(
            autoexpand=True,
            l=20,
            r=20,
            t=50, ),
        transition_duration=500,
        plot_bgcolor='white',
        title_x=0.5,
    )


def table_type(df_column):
    """
    Ensures column types are correct so filtering in main table works
    """
    # Note - this only works with Pandas >= 1.0.0

    if sys.version_info < (3, 0):  # Pandas 1.0.0 does not support Python 2
        return 'any'

    if isinstance(df_column.dtype, pd.DatetimeTZDtype):
        return 'datetime',
    elif (isinstance(df_column.dtype, pd.StringDtype) or
          isinstance(df_column.dtype, pd.BooleanDtype) or
          isinstance(df_column.dtype, pd.CategoricalDtype) or
          isinstance(df_column.dtype, pd.PeriodDtype)):
        return 'text'
    elif (isinstance(df_column.dtype, pd.SparseDtype) or
          isinstance(df_column.dtype, pd.IntervalDtype) or
          isinstance(df_column.dtype, pd.Int8Dtype) or
          isinstance(df_column.dtype, pd.Int16Dtype) or
          isinstance(df_column.dtype, pd.Int32Dtype) or
          isinstance(df_column.dtype, pd.Int64Dtype)):
        return 'numeric'
    else:
        return 'any'


def create_small_table(df_summary, name):
    """
    Args:
        df_summary (Pandas DataFrame): _description_
        name (String): unique name for DataTable id

    Returns:
        dash DataTable
    """

    table = dash_table.DataTable(
        id=f"{name}-table",
        columns=[{"name": i, "id": i} for i in df_summary.columns],
        data=df_summary.to_dict("records"),
        sort_action="native",
        sort_mode="native",
        export_format="csv",
        style_cell={
            'whiteSpace': 'normal',
            'height': 'auto',
            'padding': 5,
            'fontSize': 14,
            'font-family': 'sans-serif'
        },
    )
    return table


def create_data_table(df):
    """Create Dash datatable from Pandas DataFrame."""
    table = dash_table.DataTable(
        id="database-table",
        columns=[{'name': i, 'id': i, 'type': table_type(
            df[i])
                  } for i in df.columns],
        data=df.to_dict("records"),
        filter_action='native',
        sort_action="native",
        sort_mode="native",
        export_format="csv",
        style_header={
            'backgroundColor': '#525b63',
            'color': 'white'
        },
        style_cell={
            'whiteSpace': 'normal',
            'height': 'auto',
            'padding': 5,
            'fontSize': 14,
            'font-family': 'sans-serif'
        },
        style_table={'overflowX': 'auto'},
        page_size=10,
    )
    return table
