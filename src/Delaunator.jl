module Delaunator

using OffsetArrays

const EPSILON = 2^-52 # eps()
const EDGE_STACK = OffsetVector{UInt32}(undef,0:511)

# monotonically increases with real angle, but doesn't need expensive trigonometry
function pseudoAngle(dx, dy)
    p = dx / (abs(dx) + abs(dy))
    return (dy > 0 ? 3 - p : 1 + p) / 4 # [0..1]
end

function dist(ax, ay, bx, by)
    dx = ax - bx
    dy = ay - by
    return dx * dx + dy * dy
end

# return 2d orientation sign if we're confident in it through J. Shewchuk's error bound check
function orientIfSure(px, py, rx, ry, qx, qy)
    l = (ry - py) * (qx - px)
    r = (rx - px) * (qy - py)
    return abs(l - r) >= 3.3306690738754716e-16 * abs(l + r) ? l - r : 0
end

# a more robust orientation test that's stable in a given triangle (to fix robustness issues)
function orient(rx, ry, qx, qy, px, py)
    return (orientIfSure(px, py, rx, ry, qx, qy) ||
        orientIfSure(rx, ry, qx, qy, px, py) ||
        orientIfSure(qx, qy, px, py, rx, ry)) < 0
end

function inCircle(ax, ay, bx, by, cx, cy, px, py)
    dx = ax - px
    dy = ay - py
    ex = bx - px
    ey = by - py
    fx = cx - px
    fy = cy - py

    ap = dx * dx + dy * dy
    bp = ex * ex + ey * ey
    cp = fx * fx + fy * fy

    return dx * (ey * cp - bp * fy) -
           dy * (ex * cp - bp * fx) +
           ap * (ex * fy - ey * fx) < 0
end

function circumradius(ax, ay, bx, by, cx, cy)
    dx = bx - ax
    dy = by - ay
    ex = cx - ax
    ey = cy - ay

    bl = dx * dx + dy * dy
    cl = ex * ex + ey * ey
    d = 0.5 / (dx * ey - dy * ex)

    x = (ey * bl - dy * cl) * d
    y = (dx * cl - ex * bl) * d

    return x * x + y * y
end

function circumcenter(ax, ay, bx, by, cx, cy)
    dx = bx - ax
    dy = by - ay
    ex = cx - ax
    ey = cy - ay

    bl = dx * dx + dy * dy
    cl = ex * ex + ey * ey
    d = 0.5 / (dx * ey - dy * ex)

    x = ax + (ey * bl - dy * cl) * d
    y = ay + (dx * cl - ex * bl) * d

    return (x, y)
end

function quicksort(ids, dists, left, right)
    if right - left <= 20
        i = i = left + 1
        while i <= right
            temp = ids[i]
            tempDist = dists[temp]
            j = i - 1
            while j >= left && dists[ids[j]] > tempDist
                #ids[j + 1] = ids[j--]
            end
            ids[j + 1] = temp
            i += 1
        end
    else
        median = (left + right) >> 1
        i = left + 1
        j = right
        swap(ids, median, i)
        dists[ids[left]] > dists[ids[right]] && swap(ids, left, right)
        dists[ids[i]] > dists[ids[right]] && swap(ids, i, right)
        dists[ids[left]] > dists[ids[i]] && swap(ids, left, i)

        temp = ids[i]
        tempDist = dists[temp]
        while true
            # TODO
            #do i++; while (dists[ids[i]] < tempDist);
            #do j--; while (dists[ids[j]] > tempDist);
            j < i && break
            swap(ids, i, j)
        end
        ids[left + 1] = ids[j]
        ids[j] = temp

        if (right - i + 1 >= j - left)
            quicksort(ids, dists, i, right)
            quicksort(ids, dists, left, j - 1)
        else
            quicksort(ids, dists, left, j - 1)
            quicksort(ids, dists, i, right)
        end
    end
end

function swap(arr, i, j)
    tmp = arr[i]
    arr[i] = arr[j]
    arr[j] = tmp
end

function defaultGetX(p)
    return p[1]
end

function defaultGetY(p)
    return p[2]
end


end # module
