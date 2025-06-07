def validate_comma_separated(ctx, param, value):
    print(value)
    if not value:
        return []
    return value.split(',')