

def get_point(mat, p):
    if not point_exists(mat, p):
        return None
    return mat[p[0]][p[1]][p[2]]


def set_point(mat, p, dp):
    if not point_exists(mat, p):
        return
    mat[p[0]][p[1]][p[2]] = dp


def get_temp(mat, p):
    if not point_exists(mat, p):
        return 0
    return mat[p[0]][p[1]][p[2]].get_t()


def set_temp(mat, p, t):
    if not point_exists(mat, p):
        return
    mat[p[0]][p[1]][p[2]].set_t(t)





def calculate_dt(mat, p, dt, r):
    curr_point = get_point(mat, p)
    adj_points = get_adj(mat, p)

    curr_t = curr_point.get_t()
    cp = curr_point.get_cp()
    p = curr_point.get_p()

    a = 1 - (dt/(p*cp))*(len(adj_points)/r**3)
    b = (dt/(p*cp))*(1/r**3)
    sum_t = sum([get_temp(mat, i) for i in adj_points])
    new_t = curr_t*a + b*sum_t

    return new_t, adj_points


def get_adj(mat, p):
    ret = []
    x = p[0]
    y = p[1]
    z = p[2]

    if point_exists(mat, (x+1, y, z)):
        ret.append((x+1, y, z))
    if point_exists(mat, (x-1, y, z)):
        ret.append((x-1, y, z))
    if point_exists(mat, (x, y+1, z)):
        ret.append((x, y+1, z))
    if point_exists(mat, (x, y-1, z)):
        ret.append((x, y-1, z))
    if point_exists(mat, (x, y, z+1)):
        ret.append((x, y, z+1))
    if point_exists(mat, (x, y, z-1)):
        ret.append((x, y, z-1))
    return ret


def point_exists(mat, p):
    x = p[0]
    y = p[1]
    z = p[2]

    return not((x < 0 or x >= mat.shape[0]) or (y < 0 or y >= mat.shape[1]) or (z < 0 or z >= mat.shape[2]))



