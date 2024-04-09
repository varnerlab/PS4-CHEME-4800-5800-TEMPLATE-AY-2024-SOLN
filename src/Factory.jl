"""
    build(modeltype::Type{MySimpleProblemModel}, data::NamedTuple) -> MySimpleProblemModel 

Build a model of `MySimpleProblemModel` type using the data provided in the `data` argument (a NamedTuple).

### Arguments
- `modeltype::Type{MySimpleProblemModel}`: the type of the model to build.
- `data::NamedTuple`: the data to use to build the model.
"""
function build(modeltype::Type{MySimpleProblemModel}, data::NamedTuple)::MySimpleProblemModel
    
    # build an empty model -
    model = modeltype();
    
    # add stuff to the model from the data arg -
    model.parameters = data.parameters;
    model.initial_conditions = data.initial_conditions;
    model.time_span = data.time_span;
    model.model = data.model;
    
    # return -
    return model;
end