function quicksort(ids, dists, left, right)
    
    if right - left <= 20
        i = i = left + 1
        while i <= right
            temp = ids[i]
            tempDist = dists[temp]
            j = i - 1
            while j >= left && dists[ids[j]] > tempDist
                ids[j + 1] = ids[j]
                j -= 1
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
            i += 1
            while dists[ids[i]] < tempDist
                i += 1
            end
            j -= 1
            while dists[ids[j]] > tempDist
                j -= 1
            end
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

Base.@propagate_inbounds function swap(arr, i, j)
    tmp = arr[i]
    arr[i] = arr[j]
    arr[j] = tmp
end

function parallelquicksort(ids, dists, left, right, parallel=true)
    if right-left <= 1 <<  9 
        parallel = false 
    end 
    if right - left <= 20
        i = i = left + 1
        while i <= right
            temp = ids[i]
            tempDist = dists[temp]
            j = i - 1
            while j >= left && dists[ids[j]] > tempDist
                ids[j + 1] = ids[j]
                j -= 1
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
            i += 1
            while dists[ids[i]] < tempDist
                i += 1
            end
            j -= 1
            while dists[ids[j]] > tempDist
                j -= 1
            end
            j < i && break
            swap(ids, i, j)
        end
        ids[left + 1] = ids[j]
        ids[j] = temp

        if (right - i + 1 >= j - left)
            ii, ij, iright, ileft = i, j, right, left
            if parallel
                t = Base.Threads.@spawn quicksort(ids, dists, ii, iright, parallel)
                quicksort(ids, dists, ileft, ij - 1)
                wait(t)
            else
                quicksort(ids, dists, ii, iright)
                quicksort(ids, dists, ileft, ij - 1)
            end
        else
            ii, ij, iright, ileft = i, j, right, left
            if parallel 
                t = Base.Threads.@spawn quicksort(ids, dists, ileft, ij - 1, parallel)
                quicksort(ids, dists, ii, iright)
                wait(t)
            else 
                quicksort(ids, dists, ileft, ij - 1)
                quicksort(ids, dists, ii, iright)
            end 
        end
    end
end