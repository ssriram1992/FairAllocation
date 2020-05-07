# Note: We are using the IHME and ODE projections of April 9th 2020, so that
# they cover the same range of days.


import datetime


def state_indices():
    """ Lists the states in alphabetical order.

    Returns:
      A list of the states in alphabetical order.
    """
    with open('data/bertsimas/state_list.csv', 'r') as f:
        data = [x.strip('"\n') for x in f.readlines()[1:]]
    f.close()
    data.sort()
    return data


def ihme_clean_data():
    """ Basic clean-up of IHME data.

    Returns:
      Cleaned-up data.
    """
    with open('data/bertsimas/projected_demand_ihme.csv', 'r') as f:
        data = [[x[0].strip('"\n'), x[1], int(x[4].strip('\n'))] for x in
                [y.split(',') for y in f.readlines()[1:]]]
    f.close()
    return data


def ode_clean_data():
    """ Basic clean-up of ODE data.

    Returns:
      Cleaned-up data.
    """
    with open('data/bertsimas/projected_demand_ode.csv', 'r') as f:
        data = [[x[0].strip('"\n'), x[1], int(x[7].strip('\n'))] for x in
                [y.split(',') for y in f.readlines()[1:]]]
    f.close()
    return data


def parse_data(data):
    """ Parses cleaned IHME or ODE data.

    Arguments:
      data: Clean IHME or ODE data.
    
    Returns:
      Parsed data.
    """
    dates = list(set([x[1] for x in data]))
    dates.sort(key=lambda date: datetime.datetime.strptime(date, '%Y-%m-%d'))
    states = state_indices()
    projection = [[None for _ in range(len(states))] for _ in range(len(dates))]
    for entry in data:
        projection[dates.index(entry[1])][states.index(entry[0])] = entry[2]
    return projection


def ihme_data():
    """ Get the IHME ventilator projection.

    Returns:
      For each day, a list of ventilator demand for the states,
      indexed in alphabetical order.
    """
    return parse_data(ihme_clean_data())


def ode_data():
    """ Get the ODE ventilator projection.

    Returns:
      For each day, a list of ventilator demand for the states,
      indexed in alphabetical order.
    """
    return parse_data(ode_clean_data())


def initial_supply():
    """ Get the initial supply of ventilators.

    Returns:
      A list of the initial supply of ventilators for each state,
      indexed in alphabetical order.
    """
    with open('data/bertsimas/initial_supply.csv', 'r') as f:
        data = [[x[0].strip('"\n'), int(x[5].strip('\n'))] for x in
                [y.split(',') for y in f.readlines()[1:]]]
    f.close()
    data.sort(key=lambda entry: entry[0])
    return [x[1] for x in data]
