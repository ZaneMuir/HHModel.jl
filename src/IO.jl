function _check_filename(filename, ext)
    if splitext(filename)[2] != ext
        filename * ext
    else
        filename
    end
end

# export as atf file
function export_as_atf(filename, tspan, currents, voltages)
    _SweepStartTimesMS = "\"SweepStartTimesMS=" * join([string((idx-1)*1000.0) for idx = 1:size(currents, 2)], ",") * "\""
    _Signals = "\"Signals=\"\t" * join(["\"I_MTest 2\"\t\"I_MTest 3\"" for idx = 1:size(currents, 2)], "\t")
    _header = "\"Time (s)\"\t" * join(["\"Trace #1 (pA)\"\t\"Trace #1 (mV)\"" for idx = 1:size(currents, 2)], "\t")

    _file_header = @sprintf("ATF\t1.0\n8\t%d\n\"AcquisitionMode=Episodic Stimulation\"\n\"Comment=\"\n\"YTop=20000,1000\"\n\"YBottom=-20000,-1000\"\n\"SyncTimeUnits=5\"\n%s\n\"SignalsExported=I_MTest 2,I_MTest 3\"\n%s\n%s\n", size(currents, 2)*2+1, _SweepStartTimesMS, _Signals, _header)
    # _file_header |> println

    _format_trace = (_c, _v) -> begin
        [(@sprintf "%.5f\t%.5f\t" _c[_idx] _v[_idx]) for _idx = 1:length(_c)]
    end

    _trace_str = Vector{String}(undef, length(tspan))
    for (idx, _tstep) = enumerate(tspan)
        _data = *(_format_trace(currents[idx, :], voltages[idx, :])...)
        _trace_str[idx] = (@sprintf "%.5f\t%s\n" _tstep/1e3 _data)
    end
    _atf_str = *(_file_header, _trace_str...);

    open(_check_filename(filename, ".atf"), "w") do atffile
        write(atffile, _atf_str)
    end
    true
end

# export as MAT
function export_as_mat(filename, tspan, currents, voltages)
    MAT.matwrite(_check_filename(filename, ".mat"), Dict("meta"=>Dict(),
            "tspan" => tspan,
            "current" => currents,
            "voltage" => voltages))
    true
end

function export_as_hdf5(filename, tspan, currents, voltages; _mode="VC", _note="")
    _frequency = Int(div(1000, (tspan[2] - tspan[1])))
    h5open(_check_filename(filename, ".h5"), "w") do hfile
        _g_meta = g_create(hfile, "meta")
        _g_meta["_type"] = "SweepBasedData"
        _g_meta["_origin_format"] = "HHModel simulation"
        _g_meta["date"] = Dates.format(Dates.now(), "yyyy_mm_dd")
        _g_meta["mode"] = _mode
        _g_meta["note"] = _note

        hfile["frequency"] = _frequency
        hfile["current", "compress", 6] = currents
        hfile["voltage", "compress", 6] = voltages
    end
    true
end