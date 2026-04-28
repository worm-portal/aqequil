import pandas as pd
import copy
import math
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import aqequil

def react_water_rock(speciation, rock_composition,
                     db_next, end_T_C=300, aux_basis=[],
                     minimum_molality=1E-18, prepare_reaction_kwargs=None):
    '''
    rock_composition : dict, str, or pandas.DataFrame
        Object that contains rock composition information. For example,
        a dict with a "name" key:
        {"name": "basalt", "Si": 8.613522146, "Al": 2.956783674,
         "Fe": 1.488675742, "Mg": 1.928534280, "Ca": 2.082788925,
         "Na": 0.923200310, "K":  0.034835642, "O": 27.64123673}
        Or a pandas dataframe with columns for "name" and elements and
        rows for various rock compositions, or the filepath of a CSV
        file formatted in the same way.
    '''
    
    # check that there is only one fluid sample in the speciation that will
    # be reacted with all rock compositions.
    try:
        assert len(speciation.sample_data.keys()) == 1
    except:
        # raise exception here
        print("Could not proceed because react_water_rock() requires a speciation "
              "object containing only one fluid sample. The speciation object "
              "provided has", len(speciation.sample_data.keys()), "samples.")
        return

    sample_name = list(speciation.sample_data.keys())[0] # speciation should only have one sample
    sample_data = speciation.sample_data[sample_name]
    start_T_C = sample_data["temperature"]
    pressure_bar = sample_data["pressure"]

    ### Create a dataframe of the original fluid sample input file
    report_df = speciation.lookup(speciation.lookup('input')).copy()
    
    # Reconstruct the MultiIndex columns with cleaned level 0 names
    col_level_0 = [col[0].replace("_(input)", "", 1) for col in report_df.columns]
    col_level_0 = [col if col != "pH" else "H+" for col in col_level_0]
    col_level_1 = [col[1] for col in report_df.columns]
    
    report_df.columns = pd.MultiIndex.from_arrays([col_level_0, col_level_1])
    
    # Reset index and name it 'Sample'
    report_df.index.name = 'Sample'
    original_fluid_input_df = report_df.reset_index()

    ### Handle import of rock compositions
    if isinstance(rock_composition, str):
        rock_composition = pd.read_csv(rock_composition)

    if isinstance(rock_composition, pd.DataFrame):
        rock_comp_list = rock_composition.to_dict(orient='records')
        rock_comp_elements = [e for e in rock_composition.columns if e != "name"]
    elif isinstance(rock_composition, dict):
        rock_comp_list = [rock_composition]
        rock_comp_elements = [e for e in list(rock_composition.keys()) if e != "name"]

    ### react fluid with one or more rock compositions
    wr_speciation_list = []
    wr_input_df_list = []
    for rock_comp_dict in rock_comp_list:

        rock_name = rock_comp_dict["name"]
        rock_comp = {e:rock_comp_dict[e] for e in rock_comp_elements}

        ### define special reactant
        SR = aqequil.Reactant(
                reactant_name=rock_name,
                reactant_type="Special reactant",
                special_reactant_dict=rock_comp,
                )

        if prepare_reaction_kwargs == None:
            prepare_reaction_kwargs = {
                "t_option":1, # linear tracking with Xi
                "t_value_1":start_T_C, # starting temperature
                "t_value_2":end_T_C-start_T_C, # how much temperature should change
                "p_option":2, # constant pressure
                "p_value_1":pressure_bar, # bars pressure
                "max_n_steps":2000,
                "xi_print_int":0.1,
                "physical_system_model":"titration",
                "permit_solid_solutions":True,
                "fluid_mixing_setup":True,
                "max_n_NR_iter":99,
            }
        
        r = aqequil.Prepare_Reaction(
                reactants=[SR],
                **prepare_reaction_kwargs,
                )

        speciation_wr = copy.deepcopy(aqequil.react(speciation, r))
        speciation_wr.rename_samples({sample_name: rock_name})
        wr_speciation_list.append(speciation_wr)
        
        sp_wr_input_df = speciation_wr.create_input_file(
                db=db_next,
                filename=rock_name+".csv",
                aux_basis=aux_basis,
                minimum_molality=minimum_molality,
                return_df=True)

        wr_input_df_list.append(sp_wr_input_df)
    
    ### join original fluid and water-rock fluid dataframes into a CSV
    df_joined = aqequil.join_input_files(
          [original_fluid_input_df]+wr_input_df_list,
          fill_value="1E-18", return_df=True)

    df_joined_clean = aqequil.drop_min_molal_bases(df_joined, db=db_next, minimum_molality=1E-18)

    return df_joined_clean



def mix(speciation, fluid_1, fluid_2):

    speciation_for_mix = copy.deepcopy(speciation)
    
    # Create mixing fluid from the vent fluid result
    Mix = aqequil.Mixing_Fluid(speciation=speciation_for_mix,
                               sample_name=fluid_2,
                               mix_with=fluid_1)
    
    # Set up and run the mixing calculation
    r_mix = aqequil.Prepare_Reaction(reactants=[Mix],
                                     mineral_suppression_option="All",
                                     xi_print_int=0.1)
    
    speciation_mixed_1 = aqequil.react(speciation_for_mix, r_mix)
    m1 = speciation_mixed_1.mt(fluid_1)

    Mix = aqequil.Mixing_Fluid(speciation=speciation_for_mix,
                               sample_name=fluid_1,
                               mix_with=fluid_2)
    
    # Set up and run the mixing calculation
    r_mix = aqequil.Prepare_Reaction(reactants=[Mix],
                                     mineral_suppression_option="All",
                                     xi_print_int=0.1)
    
    speciation_mixed_2 = aqequil.react(speciation_for_mix, r_mix)
    m2 = speciation_mixed_2.mt(fluid_2)
    
    joined_mix = aqequil.join_mixes(m1, m2)

    return joined_mix



def get_axis_refs(row, col, n_cols):
    """Return the xref/yref strings for a given subplot position."""
    idx = (row - 1) * n_cols + col
    if idx == 1:
        return "x", "y"
    return f"x{idx}", f"y{idx}"
    

def combine_plots(plot_list, plot_titles=None, xlab=None, ylab=None,
                  plot_height=300, plot_width=700, title=None,
                  plot_margins=dict(l=70, r=40, t=60, b=60),
                  rows=1, cols=None, shared_yaxes=True, shared_xaxes=False,
                  plotly_kwargs=None):

    if cols==None:
        cols = math.ceil(len(plot_list)/rows)
    
    make_subplots_kwargs = {"rows":rows, "cols":cols,
                            "shared_yaxes":shared_yaxes,
                            "shared_xaxes":shared_xaxes}
    
    if plot_titles==None:
        plot_titles=[]
        for fig in plot_list:
            if fig.layout["title"]["text"] != None:
                plot_titles.append(fig.layout["title"]["text"])
            else:
                plot_titles=None
                
    if isinstance(plot_titles, list):
        make_subplots_kwargs["subplot_titles"] = tuple(plot_titles)

    if xlab==None:
        if plot_list[0].layout["xaxis"]["title"]["text"] != None:
            xlab=plot_list[0].layout["xaxis"]["title"]["text"]
        else:
            xlab=None
            
    if ylab==None:
        if plot_list[0].layout["yaxis"]["title"]["text"] != None:
            ylab=plot_list[0].layout["yaxis"]["title"]["text"]
        else:
            ylab=None
    
    if isinstance(plotly_kwargs, dict):
        plotly_kwargs = {k:v for k,v in zip(plotly_kwargs.keys(), plotly_kwargs.values()) if k not in list(make_subplots_kwargs.keys())}
        make_subplots_kwargs = make_subplots_kwargs | plotly_kwargs

    combined = make_subplots(**make_subplots_kwargs)

    colors = px.colors.qualitative.Plotly

    n_cols = make_subplots_kwargs["cols"]
    for i, fig in enumerate(plot_list):
        row, col = divmod(i, n_cols)
        row += 1  # plotly is 1-indexed
        col += 1
        xref, yref = get_axis_refs(row, col, n_cols)

        fig = copy.deepcopy(fig)
        
        for trace in fig.data:
            #trace.update(line_color=colors[i % len(colors)])
            if i > 0:
                trace.showlegend = False
            combined.add_trace(trace, row=row, col=col)
            
        # Remap and add shapes
        for shape in fig.layout.shapes:
            new_shape = shape.to_plotly_json()
    
            # Remap axis refs, preserving 'domain' suffix if present
            # e.g. 'x domain' -> 'x2 domain', 'y' -> 'y2'
            new_shape["xref"] = new_shape.get("xref", "x").replace("x", xref, 1)
            new_shape["yref"] = new_shape.get("yref", "y").replace("y", yref, 1)
    
            combined.add_shape(new_shape)
    
        combined.update_layout(height=plot_height, width=plot_width,
                              margin=plot_margins,  # left, right, top, bottom in pixels
                              template="simple_white",
                              )

    if title != None:
        # title
        combined.add_annotation(
            text=title,
            x=0.5, y=1.15,
            xref="paper", yref="paper",
            showarrow=False,
            font=dict(size=18, family="Arial"),
            xanchor="center", yanchor="bottom"
        )

    # x axis label
    combined.add_annotation(
        text=xlab,
        x=0.5, y=-0.28,
        xref="paper", yref="paper",
        showarrow=False,
        font=dict(size=14)
    )

    # y axis label
    combined.add_annotation(
        text=ylab,
        x=-0.12, y=0.5,
        xref="paper", yref="paper",
        showarrow=False,
        font=dict(size=14),
        textangle=-90
    )

    combined.show()