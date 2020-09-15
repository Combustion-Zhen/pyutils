def params2name(params, s1='_', s2='-'):
    params_str = []
    for k, v in params.items():
        try:
            param_str = '{0}{1}{2:g}'.format(k,s2,v)
        except ValueError:
            param_str = '{0}{1}{2}'.format(k,s2,v)
        params_str.append(param_str)
    return s1.join(params_str)

def name2params(name, s1='_', s2='-', default=None):
    params = {}
    params_str = name.split(s1)
    for param_str in params_str:
        param = param_str.split(s2,maxsplit=1)
        key = param[0]
        try:
            params[key] = float( param[1] )
        except ValueError:
            params[key] = param[1]
        except IndexError:
            params[key] = default
    return params

def add_suffix(name, suffix):

    if name.endswith(suffix):
        return name

    return '.'.join([name, suffix])