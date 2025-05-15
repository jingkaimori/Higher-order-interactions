function quad_clustering(inci_mat::BitMatrix)
    (n1,n2) = size(inci_mat)
    # 输出参数
    clu_coeff = zeros(Float64, n1)

    Threads.@threads for i in 1:n1
        denominator = 0
        numerator = 0
        numofedges = count(inci_mat[i, :])

        if numofedges < 2
            clu_coeff[i] = 0.0
            continue
        end

        edgelist1 = Int[]
        for a in 1:n2
            if inci_mat[i, a]
                push!(edgelist1, a)
            end
        end

        for a in 1:(numofedges - 1)
            if count(inci_mat[:, edgelist1[a]]) == 1
                continue
            end
            i1 = findfirst(x -> x, inci_mat[:, edgelist1[a]])
            i2 = findlast(x -> x, inci_mat[:, edgelist1[a]])

            for b in (a + 1):numofedges
                if count(inci_mat[:, edgelist1[b]]) == 1
                    continue
                end
                i3 = findfirst(x -> x, inci_mat[:, edgelist1[b]])
                i4 = findlast(x -> x, inci_mat[:, edgelist1[b]])

                denominator += min(count(inci_mat[:, edgelist1[a]]) - 1, count(inci_mat[:, edgelist1[b]]) - 1)

                for j in max(i1, i3):min(i2, i4)
                    if inci_mat[j, edgelist1[a]] && inci_mat[j, edgelist1[b]] && j != i
                        numerator += 1
                    end
                end
            end
        end

        if denominator == 0
            clu_coeff[i] = 0.0
        else
            clu_coeff[i] = numerator / denominator
        end
    end

    return clu_coeff
end