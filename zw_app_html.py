import dash.html as html
import dash.dcc as dcc
import dash_bootstrap_components as dbc

###############################
#   dashboard webpage setup   #
###############################
def make_dashboard_webpage(
        init_orbit_fig, init_orbit_data,
        init_pot_fig,
        init_gw_cross_fig, init_gw_plus_fig,
        init_gw_data,
        default_angmom_str,
        default_energy_str,
):
    dashboard_page = html.Div(
        children=[
            #
            #=== These apply to all rows & columns
            #
            # Interval for continuous clientside update: animation
            dcc.Interval(id="orbit-interval",
                         interval=20, # microseconds
                         disabled=False,
                         n_intervals=0,
                         max_intervals=-1),
            # two client-stored json data sets: offset and (x,y) data
            dcc.Store(id='stored-offset',
                      storage_type='memory',
                      data=0,
                      clear_data=False,
                      modified_timestamp=-1),
            dcc.Store(id='stored-orbit',
                      storage_type='memory',
                      data=0,
                      #data=init_orbit_data,                      
                      clear_data=False,
                      modified_timestamp=-1),
            dcc.Store(id='stored-energy',
                      storage_type='memory',
                      data=0,
                      clear_data=False,
                      modified_timestamp=-1),
            dcc.Store(id='stored-gw-data',
                      storage_type='memory',
                      data=0,
                      #data=init_gw_data,
                      clear_data=False,
                      modified_timestamp=-1),
            #
            #=== First Row
            #
            dbc.Row([
                # first row: (1) of (3) columns
                dbc.Col([
                    dbc.Card(
                        [
                            dbc.CardHeader("Orbit"),
                            dbc.CardBody(
                                [
                                    dcc.Graph(
                                        id='orbitgraph',
                                        figure=init_orbit_fig,
                                        style={
                                            'width': '100%',
                                            #'display': 'inline-block',
                                            #'width': '40vw'}),
                                        },
                                        config={
                                            'displayModeBar': False,
                                        }
                                    )
                                ]
                            )
                        ]
                    ),
                ], width=5),
                # first row: (2) of (3) columns
                dbc.Col([
                    dbc.Card(
                        [
                            dbc.CardHeader("Controls"),
                            dbc.CardBody(
                                [
                                    html.Button('re-calculate', id='recalculate-button',
                                                n_clicks=0,
                                                style={'display': 'inline-block'}),
                                    #'width': '10vw'}),
                                ]
                            ),
                        ]
                    )
                ], width=2),
                # first row: (3) of (3) columns
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader("Effective Potential"),
                        dbc.CardBody([
                            dbc.Row(dbc.Col(
                                html.Div(
                                    children=[
                                        dcc.Graph(
                                            id='potentialgraph',
                                            figure=init_pot_fig,
                                            style={
                                                'width': '100%',
                                                #'display': 'inline-block',
                                                #'width': '40vw'}),
                                            },
                                            config={
                                                'displayModeBar': False,
                                            }
                                        ),
                                    ]
                                )
                            )),
                            dbc.Row(dbc.Col([
                                html.Div(
                                    children=[
                                        html.Div(
                                            children="Angular Momentum: ",
                                            style={'display': 'inline-block',
                                                   'margin-right': 20}
                                        ),
                                        dcc.Input(
                                            id='angmom-val-str', type='text',
                                            value = default_angmom_str,
                                            placeholder = default_angmom_str,
                                            # don't allow update until user enters
                                            debounce = True,
                                            size = '7',
                                            style={'display': 'inline-block',
                                                   'margin-right': 40}
                                        ),
                                        html.Div(
                                            children="Energy: ",
                                            style={'display': 'inline-block',
                                                   'margin-right': 20}
                                        ),
                                        dcc.Input(
                                            id='energy-val-str', type='text',
                                            value = default_energy_str,
                                            placeholder = default_energy_str,
                                            # don't allow update until user enters
                                            debounce = True,
                                            size = '7',
                                            style={'display': 'inline-block'}
                                        )
                                    ],
                                )
                            ]))
                        ])
                    ])
                ],
                        width=5,
                )
            ],
                    style={"height": "70%"}
            ),
            #
            #=== Second Row
            #
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader("Gravitational Wave Signal"),
                        dbc.CardBody(
                            [
                                dcc.Graph(
                                    id='gw_plus_graph',
                                    figure=init_gw_plus_fig,
                                    style={
                                        #'width': '100%',
                                        'display': 'inline-block',
                                        'width': '40vw',
                                        'margin-right': 20
                                    },
                                    config={
                                        'displayModeBar': False,
                                    }
                                ),
                                dcc.Graph(
                                    id='gw_cross_graph',
                                    figure=init_gw_cross_fig,
                                    style={
                                        #'width': '100%',
                                        'display': 'inline-block',
                                        'width': '40vw',
                                    },
                                    config={
                                        'displayModeBar': False,
                                    }
                                ),
                            ]
                        )
                    ])
                ], width=12 )
            ], align='center'),
        ],
        style={"height": "100vh"},
    )
    return dashboard_page
###############################
# end dashboard webpage setup #
###############################
