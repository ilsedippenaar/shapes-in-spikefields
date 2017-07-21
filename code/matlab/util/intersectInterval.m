function intersection = intersect_interval(a, b)
    assert(a(2) >= a(1) && b(2) >= b(1));
    if a(1) > b(2) || a(2) < b(1)
      intersection = [];
    else
      intersection = [max(a(1), b(1)), min(a(2), b(2))];
    end
  end