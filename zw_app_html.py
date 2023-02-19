import dash.html as html
import dash.dcc as dcc
import dash_bootstrap_components as dbc

effpot_fontsize = '10pt'
effpot_fontsize_small = '8pt'

height_firstrow = '55vh'
height_secondrow = '35vh'
height_secondrow_card = '45vh'

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
            dcc.Store(id='stored-angmom',
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
            html.H3(
                children = "Schwarzschild geodesics for exploring \"Zoom-Whirl\" Orbits and their gravitational waves",
                style = {
                    'textAlign': 'center'
                },
            ),
            #
            #=== First Row
            #
            dbc.Row([
                dbc.CardGroup(
                    [
                        dbc.Card(
                            [
                                dbc.CardHeader("Orbit"),
                                dbc.CardBody([
                                    dcc.Graph(
                                        id='orbitgraph',
                                        figure=init_orbit_fig,
                                        style={
                                            'height': '90%',
                                        },
                                        config={
                                            'displayModeBar': False,
                                        }
                                    )
                                ]),
                            ],
                        ),
                        dbc.Card(
                            [
                                dbc.CardHeader("Effective Potential"),
                                dbc.CardBody([
                                    dcc.Graph(
                                        id='potentialgraph',
                                        figure=init_pot_fig,
                                        style={
                                            'height': '85%',
                                        },
                                        config={
                                            'displayModeBar': False,
                                        }
                                    ),
                                    html.Div(
                                        children="Angular Momentum: ",
                                        style={
                                            'display': 'inline-block',
                                            'margin-right': 10,
                                            'font-size' : effpot_fontsize,
                                        }
                                    ),
                                    dcc.Input(
                                        id='angmom-val-str', type='text',
                                        value = default_angmom_str,
                                        placeholder = default_angmom_str,
                                        # don't allow update until user enters
                                        debounce = True,
                                        size = '7',
                                        style={
                                            'display': 'inline-block',
                                            'margin-right': 5,
                                            'font-size' : effpot_fontsize,
                                        }
                                    ),
                                    html.Div(
                                        children=" M",
                                        style={
                                            'display': 'inline-block',
                                            'margin-right': 25,
                                            'font-size' : effpot_fontsize,
                                        }
                                    ),
                                    html.Div(
                                        children="Energy: ",
                                        style={
                                            'display': 'inline-block',
                                            'margin-right': 10,
                                            'font-size' : effpot_fontsize,
                                        }
                                    ),
                                    dcc.Input(
                                        id='energy-val-str', type='text',
                                        value = default_energy_str,
                                        placeholder = default_energy_str,
                                        # don't allow update until user enters
                                        debounce = True,
                                        size = '7',
                                        style={
                                            'display': 'inline-block',
                                            'font-size' : effpot_fontsize,
                                            'margin-right': 10,
                                        }
                                    ),
                                    html.Div(
                                        children="(\"Geometrized\" units, with G=c=1)",
                                        style={
                                            'display': 'inline-block',
                                            'font-size' : effpot_fontsize_small,
                                        }
                                    ),                                    
                                ]),
                            ],
                        )
                    ],
                    style = {'height' : height_firstrow}
                )     
            ]),
            #
            #=== Second Row
            #
            dbc.Row([
                dbc.Card(
                    [
                        dbc.CardHeader("Gravitational Wave Signal (for m/M = 1/100000)"),
                        dbc.CardBody(
                            [
                                dcc.Graph(
                                    id='gw_plus_graph',
                                    figure=init_gw_plus_fig,
                                    style={
                                        'display': 'inline-block',
                                        'width' : '50%',
                                        # should be able to say height=100%
                                        # but it does not work, so this must be
                                        # something smaller than the card height
                                        # specified below in "vh" (viewwindow height)
                                        'height' : height_secondrow,
                                    },
                                    config={
                                        'displayModeBar': False,
                                    }
                                ),
                                dcc.Graph(
                                    id='gw_cross_graph',
                                    figure=init_gw_cross_fig,
                                    style={
                                        'display': 'inline-block',
                                        'width' : '50%',
                                        'height' : height_secondrow,
                                    },
                                    config={
                                        'displayModeBar': False,
                                    }
                                ),
                            ]
                        )
                    ],
                    style = { 'height' : height_secondrow_card }
                )
            ])
        ],
    )
    return dashboard_page
###############################
# end dashboard webpage setup #
###############################
